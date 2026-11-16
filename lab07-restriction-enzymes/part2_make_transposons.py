import argparse
import random

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

IRS = {
    "IR1": "ATATGCTGGGT",
    "IR2": "AGAGGTCGAGGT",
    "IR3": "CCAATTCTCATT",
    "IR4": "GTGTTTCGGAAC",
    "IR5": "TTGCGTTTTAGG",
    "IR6": "ATGTCTTAGTG",
}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


class Builder:
    def __init__(self, rng):
        self.rng = rng
        self.parts = []
        self.length = 0
        self.ir_mask = []

    def random(self, n):
        self.add("".join(self.rng.choice("ACGT") for _ in range(n)), False)

    def add(self, text, is_ir):
        start = self.length
        self.parts.append(text)
        self.ir_mask.extend([is_ir] * len(text))
        self.length += len(text)
        return start

    def transposon(self, ir, core_len):
        s = self.add(ir, True)
        self.random(core_len)
        self.add(revcomp(ir), True)
        return s, self.length

    def text(self):
        return "".join(self.parts)


def build(seed):
    rng = random.Random(seed)
    b = Builder(rng)
    truth = []

    def record(name, s, e, case):
        truth.append({"name": name, "ir": IRS[name], "start": s + 1, "end": e, "length": e - s, "case": case})

    b.random(rng.randint(20, 40))
    s, e = b.transposon(IRS["IR1"], rng.randint(10, 40))
    record("IR1", s, e, "one after the other")
    b.random(rng.randint(15, 40))
    s, e = b.transposon(IRS["IR2"], rng.randint(100, 220))
    record("IR2", s, e, "one after the other (long)")
    b.random(rng.randint(15, 40))

    outer_s = b.add(IRS["IR3"], True)
    b.random(rng.randint(5, 20))
    s, e = b.transposon(IRS["IR4"], rng.randint(10, 35))
    record("IR4", s, e, "nested (inside IR3)")
    b.random(rng.randint(5, 20))
    b.add(revcomp(IRS["IR3"]), True)
    record("IR3", outer_s, b.length, "nested (contains IR4)")
    b.random(rng.randint(15, 40))

    s5 = b.add(IRS["IR5"], True)
    b.random(rng.randint(5, 15))
    s6 = b.add(IRS["IR6"], True)
    b.random(rng.randint(5, 15))
    b.add(revcomp(IRS["IR5"]), True)
    e5 = b.length
    b.random(rng.randint(5, 15))
    b.add(revcomp(IRS["IR6"]), True)
    e6 = b.length
    record("IR5", s5, e5, "overlapping (with IR6)")
    record("IR6", s6, e6, "overlapping (with IR5)")
    b.random(rng.randint(20, 40))

    seq = list(b.text())
    mask = b.ir_mask
    changed = True
    while changed:
        changed = False
        for t in truth:
            s, e = t["start"] - 1, t["end"]
            if s > 0 and e < len(seq) and seq[s - 1] == COMP[seq[e]]:
                pos = s - 1 if not mask[s - 1] else e
                if mask[pos]:
                    continue
                other = COMP[seq[e]] if pos == s - 1 else COMP[seq[s - 1]]
                seq[pos] = rng.choice([c for c in "ACGT" if c != other])
                changed = True
    truth.sort(key=lambda t: t["start"])
    return "".join(seq), truth


def main():
    p = argparse.ArgumentParser(description="Random DNA with transposons: one after the other, nested and overlapping.")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--out", default="transposons_dna.fasta")
    a = p.parse_args()

    seq, truth = build(a.seed)

    with open(a.out, "w") as f:
        f.write(f">transposons_dna seed={a.seed} length={len(seq)} transposons={len(truth)}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")
    with open("transposons_answer_key.txt", "w") as f:
        f.write("name\tstart\tend\tlength\tinverted_repeat\tcase\n")
        for t in truth:
            f.write(f"{t['name']}\t{t['start']}\t{t['end']}\t{t['length']}\t{t['ir']}\t{t['case']}\n")

    print("=" * 74)
    print("RANDOM DNA WITH TRANSPOSONS")
    print("=" * 74)
    print(f"Length: {len(seq)} bases   Transposons: {len(truth)}   Seed: {a.seed}")
    print(f"Saved:  {a.out} and transposons_answer_key.txt\n")
    print("Inverted repeats I chose:")
    for name, ir in IRS.items():
        print(f"  {name}  left {ir:<13} right {revcomp(ir)}")
    print(f"\n{'Name':<6}{'Start':>6}{'End':>6}{'Length':>8}   Case")
    print("-" * 74)
    for t in truth:
        print(f"{t['name']:<6}{t['start']:>6}{t['end']:>6}{t['length']:>8}   {t['case']}")


if __name__ == "__main__":
    main()
