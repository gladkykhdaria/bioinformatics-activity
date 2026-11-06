"""
Exercise 1 - Gel Electrophoresis Simulation
===========================================

What it does:
  Takes a DNA sequence from NCBI (1000-3000 nt), cuts 10 random fragments
  out of it (100-3000 bases each), stores them in a list, and simulates how
  they would move through an electrophoresis gel.

The logic behind it:
  1. Read the sequence
     The FASTA file from NCBI is read and all lines are joined into one
     string. Here: Ageratum yellow vein China virus, NC_021202.1, 2,739 nt.

  2. Take 10 random samples
     For each sample, pick a random length between 100 and 3000 (but not
     longer than the sequence itself), then a random start position, and
     copy that part. All 10 fragments go into one list (array).

  3. How a gel works
     DNA has a negative charge, so in an electric field it moves towards
     the + side. The gel is like a net: short pieces slip through it
     easily and travel far, long pieces get stuck more and stay near the
     start (the well).

  4. The formula
     In a real gel the distance is not linear with size. It follows the
     logarithm of the size:

         distance = a - b * log10(length)

     So 100 bp and 200 bp are far apart, but 2000 bp and 2100 bp are almost
     at the same place. The values a and b are picked so that 100 bp ends
     near the bottom of the gel and 3000 bp stays near the top.

  5. Ladder
     A real gel always has a "ladder" lane with fragments of known size
     (100, 200, 300, 500, 700, 1000, 1500, 2000, 3000 bp). We compare our
     bands with it to read their size.

  6. Draw the gel
     Each fragment becomes a band in its own lane at the position given by
     the formula. The picture is saved as gel.png.

Biology -> CS:
  Fragments = list of strings, size = len(), migration = a math function of
  log(size), the gel picture = a plot where y is the distance.

Run:
  python ex1_gel_electrophoresis.py
  python ex1_gel_electrophoresis.py ../genomes/GCF_000907355.1_ViralProj202889_genomic.fna.gz --seed 7
"""

import argparse
import gzip
import math
import random

import matplotlib.pyplot as plt

DEFAULT_GENOME = "../genomes/GCF_000907355.1_ViralProj202889_genomic.fna.gz"
LADDER = [100, 200, 300, 500, 700, 1000, 1500, 2000, 3000]
GEL_LENGTH = 100.0
MIN_BP, MAX_BP = 100, 3000


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


def take_samples(sequence, count, min_len, max_len, rng):
    max_len = min(max_len, len(sequence))
    samples = []
    for _ in range(count):
        length = rng.randint(min_len, max_len)
        start = rng.randint(0, len(sequence) - length)
        samples.append(sequence[start:start + length])
    return samples


def migration(length):
    b = (GEL_LENGTH * 0.85) / (math.log10(MAX_BP) - math.log10(MIN_BP))
    a = GEL_LENGTH * 0.92 + b * math.log10(MIN_BP)
    return a - b * math.log10(length)


def draw_gel(samples, path):
    lanes = len(samples) + 1
    fig, ax = plt.subplots(figsize=(1.0 * lanes + 1.5, 7))
    ax.set_facecolor("#0d1b2a")
    for lane in range(lanes):
        ax.add_patch(plt.Rectangle((lane + 0.15, -2.5), 0.7, 2.0, color="#24364d"))
    for size in LADDER:
        y = migration(size)
        ax.add_patch(plt.Rectangle((0.15, y - 0.5), 0.7, 1.0, color="#9fd3ff"))
        ax.text(-0.15, y, f"{size}", color="#b8c4d4", ha="right", va="center", fontsize=8)
    for i, s in enumerate(samples, start=1):
        y = migration(len(s))
        ax.add_patch(plt.Rectangle((i + 0.15, y - 0.6), 0.7, 1.2, color="#7dff9e"))
        ax.text(i + 0.5, y + 3.2, f"{len(s)} bp", color="#d8ffe3", ha="center", fontsize=7.5)
    ax.set_xlim(-0.9, lanes)
    ax.set_ylim(GEL_LENGTH + 2, -4)
    ax.set_xticks([i + 0.5 for i in range(lanes)])
    ax.set_xticklabels(["Ladder"] + [f"S{i}" for i in range(1, lanes)], fontsize=9)
    ax.set_yticks([])
    ax.set_title("Simulated agarose gel  (well at top, + pole at bottom)", fontsize=11)
    ax.text(-0.5, -1.5, "bp", color="#b8c4d4", ha="center", fontsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    return fig


def main():
    p = argparse.ArgumentParser(description="Simulate gel electrophoresis of 10 random DNA fragments.")
    p.add_argument("genome", nargs="?", default=DEFAULT_GENOME)
    p.add_argument("--samples", type=int, default=10)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--no-show", action="store_true")
    a = p.parse_args()

    name, seq = read_fasta(a.genome)
    rng = random.Random(a.seed)
    samples = take_samples(seq, a.samples, MIN_BP, MAX_BP, rng)

    print("=" * 62)
    print("GEL ELECTROPHORESIS SIMULATION")
    print("=" * 62)
    print("Sequence:", name)
    print("Length:  ", f"{len(seq):,} nt")
    print(f"Samples:  {len(samples)} random fragments, {MIN_BP}-{min(MAX_BP, len(seq))} bp")

    print("\nThe array (in the order they were taken):")
    print([len(s) for s in samples])

    print(f"\n{'Lane':<6}{'Length':>8}   {'Distance':>9}   Starts with")
    print("-" * 62)
    for i, s in enumerate(samples, start=1):
        print(f"S{i:<5}{len(s):>6} bp   {migration(len(s)):>7.1f} mm   {s[:25]}...")

    print("\nOrder on the gel, from the well (top) to the bottom:")
    order = sorted(range(len(samples)), key=lambda i: -len(samples[i]))
    print("  " + "  >  ".join(f"S{i + 1} ({len(samples[i])})" for i in order))
    print("\nLongest fragment stays at the top, shortest travels the furthest.")

    draw_gel(samples, "gel.png")
    print("\nGel picture saved as gel.png")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
