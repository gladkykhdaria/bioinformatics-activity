"""
Shotgun Genome Assembly
=======================

What it does:
  Takes a real genome, cuts 20,000 random pieces ("reads") of 100-150
  letters out of it, throws the original away, and rebuilds the genome
  only from the pieces. At the end it checks the rebuilt genome against
  the original.

The logic behind it:
  1. Sampling (the "sequencing machine")
     A window of random length (100-150) is placed at a random position and
     the letters inside are copied out. Repeat 20,000 times. This is what a
     real sequencer does: it can't read a whole genome in one go, only
     short pieces from random places.

  2. Coverage check
     Coverage = how many reads cover each letter of the genome.
     20,000 reads x ~125 letters / 2,739 letters = ~900x on average.
     If some letter has coverage 0, nothing can rebuild it.

  3. Clean up the reads
     Identical reads are merged into one. They carry no new information.

  4. Index the reads
     Each read is filed under its first k letters (k = 20), like words in a
     dictionary filed by their first letters. This lets us find matching
     reads instantly, without comparing every read with every other read
     (20,000 x 20,000 = 400 million comparisons).

  5. Greedy extension (finding "what they have in common")
     Start from one read; that is the first piece of the rebuilt genome
     (the "contig"). Look at the end of the contig and find a read whose
     beginning is exactly the same as the contig's end, for at least
     30 letters. That shared part is the overlap. Glue the read on and keep
     only its new letters. Among all matching reads, pick the one that adds
     the most new letters. Repeat until no read extends the end.

         contig: ...ACGTTAGGCATCCGA
         read:          GGCATCCGATTACG
                        |overlap||new|
         contig: ...ACGTTAGGCATCCGATTACG

  6. Extend to the left
     Same idea in the other direction: reverse every string, extend to the
     "right", reverse back.
     If some part of the genome was never read, the extension gets stuck
     at the gap. The app then starts a new piece (contig) from a read that
     isn't used yet, so you end up with several contigs instead of one.

  7. Check the result
     Compare the rebuilt genome with the original, letter by letter.

Why it works here:
  An overlap of 30 letters is only safe if that 30-letter piece appears
  once in the genome. This virus's longest repeated piece is 12 letters,
  so every overlap points to exactly one place. Genomes with long repeats
  (like the human genome) are much harder to assemble.

Biology -> CS:
  Sequencer = random sampling, overlap = string matching (suffix == prefix),
  the read index = hash table, greedy extension = greedy algorithm, the
  whole process = rebuilding a shredded document from many copies of it.

Run:
  python genome_assembly.py
  python genome_assembly.py ../genomes/GCF_000907355.1_ViralProj202889_genomic.fna.gz
  python genome_assembly.py genome.fna.gz --reads 20000 --min 100 --max 150 --seed 42
"""

import argparse
import gzip
import math
import random
import time

DEFAULT_GENOME = "../genomes/GCF_000907355.1_ViralProj202889_genomic.fna.gz"


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, parts = "", []
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    break
                name = line[1:]
            elif line:
                parts.append(line.upper())
    return name, "".join(parts)


def sample_reads(genome, n, min_len, max_len, rng):
    reads, starts = [], []
    for _ in range(n):
        length = rng.randint(min_len, max_len)
        start = rng.randint(0, len(genome) - length)
        reads.append(genome[start:start + length])
        starts.append((start, length))
    return reads, starts


def coverage(genome_len, starts):
    diff = [0] * (genome_len + 1)
    for start, length in starts:
        diff[start] += 1
        diff[start + length] -= 1
    cov, run = [], 0
    for i in range(genome_len):
        run += diff[i]
        cov.append(run)
    return cov


def build_index(reads, k):
    index = {}
    for i, r in enumerate(reads):
        index.setdefault(r[:k], []).append(i)
    return index


def extend_right(contig, reads, index, k, min_overlap, max_len, log=None):
    while True:
        best_gain, best_read, best_overlap = 0, None, 0
        start = max(0, len(contig) - max_len + 1)
        for p in range(start, len(contig) - min_overlap + 1):
            tail = contig[p:]
            for i in index.get(contig[p:p + k], ()):
                r = reads[i]
                gain = len(r) - len(tail)
                if gain > best_gain and r.startswith(tail):
                    best_gain, best_read, best_overlap = gain, r, len(tail)
        if best_read is None:
            return contig
        contig += best_read[best_overlap:]
        if log is not None:
            log.append((best_overlap, best_gain, len(contig)))


def assemble(reads, k, min_overlap, max_len):
    unique = sorted(set(reads), key=len, reverse=True)
    right_index = build_index(unique, k)
    reversed_reads = [r[::-1] for r in unique]
    left_index = build_index(reversed_reads, k)

    log, contigs, remaining = [], [], unique
    while remaining and len(contigs) < 100:
        contig = extend_right(remaining[0], unique, right_index, k, min_overlap, max_len, log if not contigs else None)
        contig = extend_right(contig[::-1], reversed_reads, left_index, k, min_overlap, max_len, log if not contigs else None)[::-1]
        contigs.append(contig)
        remaining = [r for r in remaining if r not in contig]
    contigs.sort(key=len, reverse=True)
    return contigs, len(unique), log


def compare(original, rebuilt):
    if rebuilt == original:
        return 100.0, 0
    pos = original.find(rebuilt[:30])
    if pos == -1:
        return 0.0, None
    aligned = original[pos:pos + len(rebuilt)]
    matched = sum(1 for a, b in zip(aligned, rebuilt) if a == b)
    return matched / len(original) * 100, len(aligned) - matched


def main():
    p = argparse.ArgumentParser(description="Simulate shotgun sequencing and rebuild a genome from reads.")
    p.add_argument("genome", nargs="?", default=DEFAULT_GENOME)
    p.add_argument("--reads", type=int, default=20000)
    p.add_argument("--min", dest="min_len", type=int, default=100)
    p.add_argument("--max", dest="max_len", type=int, default=150)
    p.add_argument("--overlap", type=int, default=30)
    p.add_argument("--k", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    a = p.parse_args()

    name, genome = read_fasta(a.genome)
    g = len(genome)
    rng = random.Random(a.seed)

    print("=" * 64)
    print("SHOTGUN GENOME ASSEMBLY")
    print("=" * 64)
    print("Genome:  ", name)
    print("Length:  ", f"{g:,} nt")

    t0 = time.time()
    reads, starts = sample_reads(genome, a.reads, a.min_len, a.max_len, rng)
    cov = coverage(g, starts)
    avg_len = sum(len(r) for r in reads) / len(reads)
    expected = a.reads * avg_len / g
    zero = sum(1 for c in cov if c == 0)

    print("\n1. SAMPLING")
    print(f"   Reads taken:          {a.reads:,} (length {a.min_len}-{a.max_len}, avg {avg_len:.1f})")
    print(f"   Expected coverage:    {expected:.0f}x  (reads x avg length / genome length)")
    print(f"   Real coverage:        min {min(cov)}x, max {max(cov)}x, mean {sum(cov) / g:.0f}x")
    print(f"   Letters never read:   {zero}")
    print(f"   First read:           {reads[0][:50]}...")

    contigs, unique_count, log = assemble(reads, a.k, a.overlap, a.max_len)
    rebuilt = contigs[0]
    t1 = time.time()

    print("\n2. ASSEMBLY")
    print(f"   Unique reads:         {unique_count:,} (duplicates removed)")
    print(f"   Contigs built:        {len(contigs)}" + ("  (1 = one complete piece)" if len(contigs) == 1 else "  (gaps: some regions were never read)"))
    if len(contigs) > 1:
        print("   Contig lengths:       " + ", ".join(f"{len(c):,}" for c in contigs[:10]) + (" ..." if len(contigs) > 10 else ""))
    print(f"   Extension steps:      {len(log)} (first contig)")
    print("   First 5 steps (overlap -> new letters -> contig length):")
    for ov, gain, length in log[:5]:
        print(f"      overlap {ov:>3} nt  +{gain:>3} nt  ->  {length:>5,} nt")
    print(f"   Time:                 {t1 - t0:.2f} s")

    identity, mismatches = compare(genome, rebuilt)
    print("\n3. RESULT")
    print(f"   Original length:      {g:,} nt")
    print(f"   Rebuilt length:       {len(rebuilt):,} nt" + (" (longest contig)" if len(contigs) > 1 else ""))
    if rebuilt == genome:
        print("   Match:                100% identical. The genome was rebuilt perfectly.")
    else:
        print(f"   Match:                {identity:.2f}% of the original, {mismatches} mismatches")
    print(f"   Rebuilt start:        {rebuilt[:50]}...")
    print(f"   Original start:       {genome[:50]}...")

    theory = g * math.exp(-a.reads * avg_len / g)
    print("\n4. THEORY (Lander-Waterman)")
    print(f"   Expected letters left uncovered: G x e^(-coverage) = {theory:.2e}")
    for n in (50, 100, 200, 500, 1000):
        c = n * avg_len / g
        print(f"   With {n:>5} reads -> coverage {c:5.1f}x -> ~{g * math.exp(-c):7.1f} letters missed")


if __name__ == "__main__":
    main()
