"""
Exercise 1 - Artificial DNA with Transposable Elements
======================================================

What it does:
  Builds a random DNA sequence of 200-400 bases and inserts 3-4
  transposable elements (transposons) into it, the same way it happens in
  a real cell. Saves the sequence as a FASTA file and the true positions
  in a separate answer key, so the detector in Exercise 2 can be checked
  without seeing the answer.

What a transposon looks like:
  A transposon is a piece of DNA that can cut itself out and jump to
  another place. Its two ends are inverted repeats (IR): the right end is
  the left end read backwards on the other strand (reverse complement).

      left IR      core            right IR
      ATGCAGTC  ...anything...  GACTGCAT

      reverse complement of ATGCAGTC:
      reverse  -> CTGACGTA
      complement each letter (A<->T, C<->G) -> GACTGCAT

How it is inserted (the picture from the board):
  1. The enzyme (transposase) cuts the host DNA at a target site, but the
     two strands are cut a few bases apart (a staggered cut).
  2. The transposon goes into the cut.
  3. The cell fills in the single-stranded gaps on both sides.
  4. Result: the short target site now appears twice, once on each side.
     These copies are the direct repeats, also called the target site
     duplication (TSD).

      host:    ...TAG ATGCA CAG...
      after:   ...TAG ATGCA [IR ... transposon ... IR] ATGCA CAG...
                      ^^^^^                             ^^^^^
                      direct repeat                     direct repeat

The logic behind the code:
  1. Pick how many transposons (3 or 4) and build each one:
     a random left IR (8-10 bases), a random core (15-35 bases), and the
     right IR = reverse complement of the left IR.
  2. Pick a total length between 200 and 400 and make random host DNA
     that leaves exactly enough space for the transposons.
  3. Split the host into 3-4 equal sections and pick one random insertion
     point inside each section, so the transposons never overlap.
  4. Insert from left to right: copy the host up to the target site,
     add the 5-base target site, the transposon, and the target site again.
  5. Write artificial_dna.fasta and answer_key.txt.

Run:
  python ex1_make_sequence.py
  python ex1_make_sequence.py --seed 7 --count 3
"""

import argparse
import random

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


def random_dna(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_transposon(rng, ir_min, ir_max, core_min, core_max):
    ir = random_dna(rng.randint(ir_min, ir_max), rng)
    while ir == revcomp(ir):
        ir = random_dna(len(ir), rng)
    core = random_dna(rng.randint(core_min, core_max), rng)
    return ir, core, ir + core + revcomp(ir)


def build(seed, count, tsd_len, total_min, total_max):
    rng = random.Random(seed)
    tes = [make_transposon(rng, 8, 10, 15, 35) for _ in range(count)]
    te_space = sum(len(te) + tsd_len for _, _, te in tes)
    gap = 10
    min_host = count * (tsd_len + 2 * gap)
    while te_space + min_host > total_max:
        tes = [make_transposon(rng, 8, 10, 15, 35) for _ in range(count)]
        te_space = sum(len(te) + tsd_len for _, _, te in tes)
    total = rng.randint(max(total_min, te_space + min_host), total_max)
    host_len = total - te_space
    host = random_dna(host_len, rng)

    section = host_len // count
    points = [rng.randint(k * section + gap, (k + 1) * section - tsd_len - gap) for k in range(count)]

    order = list(range(count))
    rng.shuffle(order)
    parts, truth, pos, last = [], [], 0, 0
    for p, ti in zip(points, order):
        ir, core, te = tes[ti]
        chunk = host[last:p + tsd_len]
        parts.append(chunk)
        pos += len(chunk)
        tsd = host[p:p + tsd_len]
        truth.append({"start": pos + 1, "end": pos + len(te), "ir": ir, "tsd": tsd, "length": len(te)})
        parts.append(te)
        pos += len(te)
        last = p
    parts.append(host[last:])
    seq = "".join(parts)
    return seq, truth


def main():
    p = argparse.ArgumentParser(description="Make an artificial DNA sequence with transposons.")
    p.add_argument("--seed", type=int, default=6)
    p.add_argument("--count", type=int, choices=[3, 4], default=4)
    p.add_argument("--tsd", type=int, default=5)
    p.add_argument("--out", default="artificial_dna.fasta")
    a = p.parse_args()

    seq, truth = build(a.seed, a.count, a.tsd, 200, 400)

    with open(a.out, "w") as f:
        f.write(f">artificial_dna seed={a.seed} length={len(seq)} transposons={a.count}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")

    with open("answer_key.txt", "w") as f:
        f.write("number\tstart\tend\tlength\tinverted_repeat\tdirect_repeat\n")
        for n, t in enumerate(truth, 1):
            f.write(f"{n}\t{t['start']}\t{t['end']}\t{t['length']}\t{t['ir']}\t{t['tsd']}\n")

    print("=" * 64)
    print("ARTIFICIAL DNA WITH TRANSPOSONS")
    print("=" * 64)
    print(f"Length:      {len(seq)} bases")
    print(f"Transposons: {a.count}")
    print(f"Saved:       {a.out} (sequence) and answer_key.txt (true positions)\n")
    print(f"{'#':<3}{'Start':>6}{'End':>6}{'Len':>5}   {'Left IR':<11}{'Right IR':<11}Direct repeat")
    print("-" * 64)
    for n, t in enumerate(truth, 1):
        print(f"{n:<3}{t['start']:>6}{t['end']:>6}{t['length']:>5}   {t['ir']:<11}{revcomp(t['ir']):<11}{t['tsd']}")

    print("\nSequence (transposons in lower case):")
    marked = list(seq)
    for t in truth:
        for i in range(t["start"] - 1, t["end"]):
            marked[i] = marked[i].lower()
    marked = "".join(marked)
    for i in range(0, len(marked), 60):
        print(f"{i + 1:>4}  {marked[i:i + 60]}")


if __name__ == "__main__":
    main()
