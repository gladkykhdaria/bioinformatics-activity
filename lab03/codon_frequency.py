"""
Codon Frequency Analyzer
========================

What it does:
  Opens a FASTA file with a viral or bacterial genome, counts how often each
  of the 64 possible codons appears, and shows the 10 most frequent ones.

The logic behind it:
  1. Read the file
     A FASTA file has a header line starting with ">" (the genome's name),
     then the DNA sequence split over many lines. We keep the name, drop the
     header, and glue all sequence lines into one long string.
     NCBI files come compressed (.fna.gz); those are unzipped on the fly
     while reading, so you can pass the .gz file as it is.

  2. Prepare 64 counters
     4 bases (A, C, G, T) in 3 positions = 4^3 = 64 codons. We build all of
     them first and set each counter to 0, so even codons that never appear
     still show up with 0.

  3. Slide a window of exactly 3 letters
     Put a 3-letter window at position 0, read the codon inside, add 1 to its
     counter. Move the window 1 letter to the right and repeat until the
     window reaches the end of the genome.

         ATGCAT...
         [ATG]        -> ATG + 1
          [TGC]       -> TGC + 1
           [GCA]      -> GCA + 1
            [CAT]     -> CAT + 1

     A genome of length N gives N - 2 windows. A window with any letter that
     is not A, C, G or T (for example N = unknown base) is skipped.

  4. Turn counts into frequencies
     frequency = count of the codon / total number of windows counted.
     All 64 frequencies together add up to 100%.

  5. Sort and show the top 10
     Sort the codons from highest to lowest count and print the first 10,
     with the amino acid each one codes for and a small text bar chart.

Biology -> CS:
  The genome is the file on the hard drive, the window is the read head,
  the 64 counters are a dictionary (hash table), and the top 10 is a sort.
  Different organisms prefer different codons for the same amino acid
  ("codon usage bias"), so this table works like a fingerprint of the genome.

Run:
  python codon_frequency.py genome.fasta
  python codon_frequency.py GCF_000907355.1_ViralProj202889_genomic.fna.gz
  python codon_frequency.py            (it will ask for the file name)
"""

import gzip
import sys
from itertools import product

genetic_code = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp",

    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",

    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",

    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}


def read_fasta(path):
    name = ""
    lines = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if not name:
                    name = line[1:]
            else:
                lines.append(line.upper())
    sequence = "".join(lines).replace("U", "T")
    return name, sequence


def count_codons(sequence):
    counts = {"".join(c): 0 for c in product("ACGT", repeat=3)}
    skipped = 0
    for i in range(len(sequence) - 2):
        window = sequence[i:i + 3]
        if window in counts:
            counts[window] += 1
        else:
            skipped += 1
    return counts, skipped


def show_top(counts, total, top=10):
    ranked = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    max_count = ranked[0][1] if ranked[0][1] else 1

    print(f"\nTop {top} codons")
    print(f"{'Rank':<6}{'Codon':<7}{'AA':<6}{'Count':>7}{'Freq %':>9}   Chart")
    print("-" * 60)
    for rank, (codon, count) in enumerate(ranked[:top], 1):
        freq = count / total * 100 if total else 0
        bar = "#" * round(count / max_count * 25)
        print(f"{rank:<6}{codon:<7}{genetic_code[codon]:<6}{count:>7}{freq:>8.2f}%   {bar}")
    return ranked


def show_all(ranked, total):
    print("\nAll 64 codons")
    print("-" * 60)
    row = []
    for codon, count in sorted(ranked):
        freq = count / total * 100 if total else 0
        row.append(f"{codon} {count:>5} ({freq:4.2f}%)")
        if len(row) == 4:
            print("   ".join(row))
            row = []


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else input("FASTA file: ").strip().strip('"')

    try:
        name, sequence = read_fasta(path)
    except OSError as e:
        print("Could not open file:", e)
        return

    if len(sequence) < 3:
        print("The sequence is too short to contain a codon.")
        return

    counts, skipped = count_codons(sequence)
    total = sum(counts.values())

    print("=" * 60)
    print("CODON FREQUENCY ANALYSIS")
    print("=" * 60)
    print("Genome:        ", name or "(no name)")
    print("Length:        ", f"{len(sequence):,} nt")
    print("Windows read:  ", f"{total:,}")
    if skipped:
        print("Windows skipped (non-ACGT):", skipped)
    gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    print("GC content:    ", f"{gc:.2f}%")

    ranked = show_top(counts, total)

    answer = input("\nShow all 64 codons? y/N: ").strip().lower()
    if answer == "y":
        show_all(ranked, total)


if __name__ == "__main__":
    main()
