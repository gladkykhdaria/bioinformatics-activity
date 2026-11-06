"""
Exercise 2 - Detecting Repeats in DNA
=====================================

What it does:
  Takes a DNA sequence from NCBI (1000-3000 nt) and finds places where a
  short unit of 3, 4 or 5 bases is repeated right after itself, so that
  the whole repeated region is 6 to 10 bases long.

      ATG ATG        unit ATG, 2 copies, 6 bases
      ATG ATG ATG    unit ATG, 3 copies, 9 bases
      GATC GATC      unit GATC, 2 copies, 8 bases
      ACGTT ACGTT    unit ACGTT, 2 copies, 10 bases

  These are called tandem repeats (or microsatellites when the unit is short).

The logic behind it:
  1. Read the sequence from the FASTA file into one string.

  2. For each unit length k = 3, 4, 5:
     Slide along the sequence and compare every letter with the letter k
     positions later. If they are the same, the pattern continues.
     A row of matches means the sequence is repeating with period k.

         position:  0 1 2 3 4 5
         letter:    A T G A T G
         compare s[i] with s[i+3]:  A=A, T=T, G=G  -> 3 matches in a row
         region length = matches + k = 3 + 3 = 6  -> ATG x 2

  3. Count full copies: copies = region length // k. Keep it only if there
     are at least 2 copies and the repeated part is 6-10 bases long.

  4. Skip "fake" units that are themselves made of a smaller repeat, like
     AAA (that is just A A A) or ATAT (that is AT AT). They belong to a
     smaller unit, not a real 3-5 base unit.

  5. Longer repeats (more than 10 bases) are listed separately so they
     don't get lost.

  6. Print all repeats with their position, unit and copies, and count
     which units appear most often.

  7. Compare with chance
     Short repeats like ATG ATG also appear in random DNA just by luck.
     To check if the virus has more repeats than expected, the same
     letters are shuffled into a random order 100 times and the repeats
     are counted again. If the real count is close to the shuffled
     average, the repeats are mostly there by chance.

Why this way:
  Comparing s[i] with s[i+k] finds each repeated region once, as one block.
  A simple "search every 3-letter word twice in a row" would find the same
  region several times (ATGATGA contains ATG ATG and also TGA TGA).

Biology -> CS:
  Tandem repeat = periodic substring, the scan = one pass per unit length,
  O(n) each, so the whole thing is very fast even for big genomes.

Run:
  python ex2_tandem_repeats.py
  python ex2_tandem_repeats.py ../genomes/GCF_000907355.1_ViralProj202889_genomic.fna.gz
"""

import argparse
import gzip
import random
from collections import Counter

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


def is_primitive(unit):
    n = len(unit)
    for d in range(1, n):
        if n % d == 0 and unit[:d] * (n // d) == unit:
            return False
    return True


def find_repeats(seq, unit_min=3, unit_max=5):
    found = []
    for k in range(unit_min, unit_max + 1):
        i = 0
        while i < len(seq) - k:
            if seq[i] != seq[i + k]:
                i += 1
                continue
            j = i
            while j < len(seq) - k and seq[j] == seq[j + k]:
                j += 1
            region = (j - i) + k
            copies = region // k
            unit = seq[i:i + k]
            if copies >= 2 and is_primitive(unit):
                found.append({
                    "start": i + 1,
                    "unit": unit,
                    "k": k,
                    "copies": copies,
                    "length": copies * k,
                    "text": seq[i:i + copies * k],
                })
            i = j + 1
    found.sort(key=lambda r: r["start"])
    return found


def main():
    p = argparse.ArgumentParser(description="Find tandem repeats with a 3-5 base unit.")
    p.add_argument("genome", nargs="?", default=DEFAULT_GENOME)
    p.add_argument("--min", dest="min_total", type=int, default=6)
    p.add_argument("--max", dest="max_total", type=int, default=10)
    a = p.parse_args()

    name, seq = read_fasta(a.genome)
    repeats = find_repeats(seq)
    inside = [r for r in repeats if a.min_total <= r["length"] <= a.max_total]
    longer = [r for r in repeats if r["length"] > a.max_total]

    print("=" * 66)
    print("TANDEM REPEAT DETECTION")
    print("=" * 66)
    print("Sequence:", name)
    print("Length:  ", f"{len(seq):,} nt")
    print(f"Looking for: unit of 3-5 bases, repeated region {a.min_total}-{a.max_total} bases")

    print(f"\nFound {len(inside)} repeats\n")
    print(f"{'#':<4}{'Position':>9}  {'Unit':<7}{'Copies':>7}{'Length':>8}   Sequence")
    print("-" * 66)
    for n, r in enumerate(inside, 1):
        shown = " ".join(r["text"][x:x + r["k"]] for x in range(0, r["length"], r["k"]))
        print(f"{n:<4}{r['start']:>9}  {r['unit']:<7}{r['copies']:>7}{r['length']:>8}   {shown}")

    by_k = Counter(r["k"] for r in inside)
    print("\nBy unit length:")
    for k in (3, 4, 5):
        print(f"  {k} bases: {by_k.get(k, 0)}")

    units = Counter(r["unit"] for r in inside)
    print("\nMost common units:")
    for unit, c in units.most_common(5):
        print(f"  {unit:<6} {c} times")

    covered = set()
    for r in inside:
        covered.update(range(r["start"], r["start"] + r["length"]))
    print(f"\nBases inside repeats: {len(covered)} of {len(seq):,} ({len(covered) / len(seq) * 100:.1f}%)")

    rng = random.Random(1)
    letters = list(seq)
    shuffled_counts = []
    for _ in range(100):
        rng.shuffle(letters)
        rnd = find_repeats("".join(letters))
        shuffled_counts.append(sum(1 for r in rnd if a.min_total <= r["length"] <= a.max_total))
    avg = sum(shuffled_counts) / len(shuffled_counts)
    print(f"\nCompared with chance (same letters shuffled, 100 times):")
    print(f"  Real sequence:     {len(inside)} repeats")
    print(f"  Shuffled average:  {avg:.1f} repeats (min {min(shuffled_counts)}, max {max(shuffled_counts)})")
    if min(shuffled_counts) <= len(inside) <= max(shuffled_counts):
        print("  -> Within the random range: most of these short repeats can appear by chance.")
    else:
        print("  -> Outside the random range: this sequence has unusually " + ("many" if len(inside) > avg else "few") + " repeats.")

    if longer:
        print(f"\nLonger than {a.max_total} bases (not counted above):")
        for r in longer:
            print(f"  position {r['start']}: {r['unit']} x {r['copies']} = {r['length']} bases")
    else:
        print(f"\nNo repeats longer than {a.max_total} bases in this sequence.")


if __name__ == "__main__":
    main()
