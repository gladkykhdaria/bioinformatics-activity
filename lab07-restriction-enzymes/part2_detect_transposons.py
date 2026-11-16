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


def build_index(seq, k):
    index = {}
    for j in range(len(seq) - k + 1):
        index.setdefault(seq[j:j + k], []).append(j)
    return index


def grow(seq, s, e, min_gap):
    n = len(seq)
    while s > 0 and e < n and seq[s - 1] == COMP[seq[e]]:
        s -= 1
        e += 1
    m = 0
    while (e - m - 1) - (s + m + 1) >= min_gap and seq[s + m] == COMP[seq[e - 1 - m]]:
        m += 1
    return s, e, m


def detect(seq, min_ir=11, min_gap=10):
    index = build_index(seq, min_ir)
    found = {}
    for i in range(len(seq) - min_ir + 1):
        for j in index.get(revcomp(seq[i:i + min_ir]), ()):
            if j - (i + min_ir) < min_gap:
                continue
            s, e, m = grow(seq, i, j + min_ir, min_gap)
            if m >= min_ir and (s, e) not in found:
                found[(s, e)] = {"start": s + 1, "end": e, "length": e - s,
                                 "ir": seq[s:s + m], "ir_right": seq[e - m:e]}
    return sorted(found.values(), key=lambda t: (t["start"], -t["end"]))


def relation(a, b):
    if a["start"] < b["start"] and b["end"] < a["end"]:
        return "contains"
    if b["start"] < a["start"] and a["end"] < b["end"]:
        return "inside"
    if a["start"] < b["start"] <= a["end"] < b["end"] or b["start"] < a["start"] <= b["end"] < a["end"]:
        return "overlaps"
    return None


def describe(found):
    for k, a in enumerate(found):
        notes = []
        for j, b in enumerate(found):
            if j != k:
                r = relation(a, b)
                if r:
                    notes.append(f"{r} #{j + 1}")
        a["relation"] = ", ".join(notes) if notes else "on its own"
    return found


def draw(found, length, width=72):
    levels = []
    rows = []
    for k, t in enumerate(found):
        a = int((t["start"] - 1) / length * width)
        b = max(a + 2, int(t["end"] / length * width))
        lvl = next((i for i, end in enumerate(levels) if end < a), None)
        if lvl is None:
            lvl = len(levels)
            levels.append(0)
        levels[lvl] = b
        rows.append((lvl, a, b, k + 1))
    lines = [[" "] * (width + 4) for _ in levels]
    for lvl, a, b, num in rows:
        lines[lvl][a] = "["
        lines[lvl][b - 1] = "]"
        for x in range(a + 1, b - 1):
            lines[lvl][x] = "-"
        label = str(num)
        mid = (a + b) // 2
        for c, ch in enumerate(label):
            if a < mid + c < b - 1:
                lines[lvl][mid + c] = ch
    print("  " + "=" * width + "  DNA")
    for row in lines:
        print("  " + "".join(row).rstrip())


def read_key(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = [l.rstrip("\n").split("\t") for l in f][1:]
    return [{"name": r[0], "start": int(r[1]), "end": int(r[2]), "case": r[5]} for r in rows if len(r) >= 6]


def main():
    p = argparse.ArgumentParser(description="Find transposons from their inverted repeats, including nested and overlapping ones.")
    p.add_argument("fasta", nargs="?", default="transposons_dna.fasta")
    p.add_argument("--min-ir", type=int, default=11)
    p.add_argument("--min-gap", type=int, default=10)
    p.add_argument("--key", default="transposons_answer_key.txt")
    a = p.parse_args()

    name, seq = read_fasta(a.fasta)
    found = describe(detect(seq, a.min_ir, a.min_gap))

    print("=" * 78)
    print("TRANSPOSON DETECTION")
    print("=" * 78)
    print("Sequence:", name)
    print("Length:  ", len(seq), "bases")
    print(f"Rule:     inverted repeats of at least {a.min_ir} bases, at least {a.min_gap} bases apart, no maximum\n")
    print(f"Found {len(found)} transposons:\n")
    print(f"{'#':<3}{'Start':>6}{'End':>6}{'Length':>8}   {'Left IR':<14}{'Right IR':<14}Position")
    print("-" * 78)
    for k, t in enumerate(found, 1):
        print(f"{k:<3}{t['start']:>6}{t['end']:>6}{t['length']:>8}   {t['ir']:<14}{t['ir_right']:<14}{t['relation']}")

    print("\nMap (nested = one bracket inside another, overlapping = brackets that cross):\n")
    draw(found, len(seq))

    key = read_key(a.key)
    if key is not None:
        got = {(t["start"], t["end"]) for t in found}
        print(f"\nCheck against {a.key}:")
        for t in key:
            ok = (t["start"], t["end"]) in got
            print(f"  {t['name']}  {t['start']:>4}-{t['end']:<4} {t['case']:<28} " + ("found exactly" if ok else "MISSED"))
        extra = got - {(t["start"], t["end"]) for t in key}
        hits = sum(1 for t in key if (t["start"], t["end"]) in got)
        print(f"  Result: {hits} of {len(key)} found exactly, {len(extra)} false hits")


if __name__ == "__main__":
    main()
