"""
Exercise 2 - Detecting Transposable Elements
============================================

What it does:
  Reads a DNA sequence from a FASTA file and finds the transposons in it,
  printing the start and end of each one. It only uses the sequence. If
  answer_key.txt is next to it, the results are checked against it at the
  end.

What it looks for (two signals together):
  1. Inverted repeats (IR): a piece of DNA at the start of the transposon,
     and its reverse complement at the end.

         ATGCAGTC ........ GACTGCAT
         |  IR  |          | IR'  |

  2. Direct repeats (target site duplication): the same few bases right
     before the start and right after the end.

         CCTAG [ATGCAGTC ...... GACTGCAT] CCTAG

  One of these alone can happen by chance in random DNA. Both together,
  at exactly the right places, almost never happen by chance.

The logic behind the code:
  1. Seed: first every 7-letter word is saved in a dictionary with its
     positions. Then for every word, its reverse complement is looked up
     in the dictionary, and only matches 24 to 150 letters to the right
     are kept. Each match is a possible pair of inverted repeats.
     (The dictionary avoids comparing every position with every other
     position.)

  2. Grow: extend the match outwards and inwards as long as the letters
     are still reverse complements. This gives the full IR length.

  3. Find the exact edges: start from the outer edge of the IR and check
     if there is a direct repeat (4-10 letters) just outside it. If not,
     move both edges one letter inwards and check again. The first edge
     with a direct repeat is the real start and end.
     This matters because sometimes the letters next to the IR match by
     chance and make the IR look one letter longer than it is. Then the
     direct repeat check fails at the wrong edge, and the next try finds
     the right one.

  4. Score: score = IR length + 2 x direct repeat length. Keep only
     candidates with IR >= 7 and direct repeat >= 4.

  5. Remove overlaps: if two candidates overlap, keep the one with the
     higher score.

  6. Print start and end (1-based, like in biology), the IR and the direct
     repeat. If there is an answer key, compare.

How likely is a false hit by chance?
  A 7-letter word has 4^7 = 16,384 possible versions. In 400 bases with a
  window of about 125, we expect about 400 x 125 / 16,384 = 3 random IR
  pairs. A 4-letter direct repeat on top of that has a chance of
  1 / 4^4 = 1/256. So we expect about 3 / 256 = 0.01 false transposons.

Run:
  python ex2_detect_transposons.py
  python ex2_detect_transposons.py my_sequence.fasta --min-ir 7 --min-tsd 4
"""

import argparse
import os

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def revcomp(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def read_fasta(path):
    with open(path, encoding="utf-8") as f:
        lines = [l.strip() for l in f if l.strip()]
    name = lines[0][1:] if lines and lines[0].startswith(">") else ""
    seq = "".join(l for l in lines if not l.startswith(">")).upper()
    return name, "".join(c for c in seq if c in "ACGT")


def tsd_length(seq, s, e, t_min, t_max):
    for t in range(t_max, t_min - 1, -1):
        if s - t >= 0 and e + t <= len(seq) and seq[s - t:s] == seq[e:e + t]:
            return t
    return 0


def detect(seq, min_ir=7, min_tsd=4, max_tsd=10, min_len=24, max_len=150):
    n = len(seq)
    comp = [COMP[c] for c in seq]
    seen, candidates = set(), []
    index = {}
    for j in range(n - min_ir + 1):
        index.setdefault(seq[j:j + min_ir], []).append(j)

    for i in range(n - min_ir + 1):
        lo = i + min_len - min_ir
        hi = min(n - min_ir, i + max_len - min_ir)
        for j in index.get(revcomp(seq[i:i + min_ir]), ()):
            if j < lo or j > hi:
                continue
            s, e = i, j + min_ir
            while s > 0 and e < n and seq[s - 1] == comp[e]:
                s -= 1
                e += 1
            m = 0
            while s + m < e - 1 - m and seq[s + m] == comp[e - 1 - m]:
                m += 1
            if (s, e) in seen:
                continue
            seen.add((s, e))

            best = None
            for d in range(0, m - min_ir + 1):
                t = tsd_length(seq, s + d, e - d, min_tsd, max_tsd)
                if t:
                    best = (s + d, e - d, m - d, t)
                    break
            if best:
                bs, be, ir, t = best
                candidates.append({
                    "start": bs + 1, "end": be, "length": be - bs,
                    "ir": seq[bs:bs + ir], "ir_right": seq[be - ir:be],
                    "tsd": seq[bs - t:bs], "score": ir + 2 * t,
                })

    candidates.sort(key=lambda c: -c["score"])
    chosen = []
    for c in candidates:
        if all(c["end"] < o["start"] or c["start"] > o["end"] for o in chosen):
            chosen.append(c)
    return sorted(chosen, key=lambda c: c["start"])


def read_key(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = [l.rstrip("\n").split("\t") for l in f][1:]
    return [(int(r[1]), int(r[2])) for r in rows if len(r) >= 3]


def main():
    p = argparse.ArgumentParser(description="Find transposons (start, end) in a DNA sequence.")
    p.add_argument("fasta", nargs="?", default="artificial_dna.fasta")
    p.add_argument("--min-ir", type=int, default=7)
    p.add_argument("--min-tsd", type=int, default=4)
    p.add_argument("--key", default="answer_key.txt")
    a = p.parse_args()

    name, seq = read_fasta(a.fasta)
    found = detect(seq, a.min_ir, a.min_tsd)

    print("=" * 70)
    print("TRANSPOSON DETECTION")
    print("=" * 70)
    print("Sequence:", name)
    print("Length:  ", len(seq), "bases")
    print(f"Rules:    inverted repeats >= {a.min_ir} bases, direct repeats >= {a.min_tsd} bases\n")
    print(f"Found {len(found)} transposable elements:\n")
    print(f"{'#':<3}{'Start':>6}{'End':>6}{'Len':>5}   {'Left IR':<12}{'Right IR':<12}Direct repeat")
    print("-" * 70)
    for k, c in enumerate(found, 1):
        print(f"{k:<3}{c['start']:>6}{c['end']:>6}{c['length']:>5}   {c['ir']:<12}{c['ir_right']:<12}{c['tsd']}")

    print("\nHow each one looks in the sequence:")
    for k, c in enumerate(found, 1):
        s, e, t = c["start"] - 1, c["end"], len(c["tsd"])
        ir = len(c["ir"])
        print(f"  {k}. {seq[s - t:s]} [{seq[s:s + ir]} ... {seq[e - ir:e]}] {seq[e:e + t]}")

    key = read_key(a.key)
    if key is not None:
        got = {(c["start"], c["end"]) for c in found}
        exact = sum(1 for k in key if k in got)
        print(f"\nCheck against {a.key}:")
        for st, en in key:
            print(f"  true {st:>4}-{en:<4}  ->  " + ("found exactly" if (st, en) in got else "MISSED"))
        extra = [g for g in got if g not in set(key)]
        print(f"  Result: {exact} of {len(key)} found with exact start and end, {len(extra)} false hits")


if __name__ == "__main__":
    main()
