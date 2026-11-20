import argparse
import os

import matplotlib.pyplot as plt

S = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
EXPECTED_CG = 29.27
EXPECTED_IC = 27.53


def cg_content(seq):
    return (seq.count("C") + seq.count("G")) / len(seq) * 100


def f(a, b):
    return 1 if a == b else 0


def self_match(seq):
    n = len(seq)
    total = 0
    for u in range(1, n):
        same = sum(f(seq[i], seq[i + u]) for i in range(n - u))
        total += same / (n - u) * 100
    return total / (n - 1)


def information_content(seq):
    return 100 - self_match(seq)


def sliding_windows(seq, size):
    return [(i, seq[i:i + size]) for i in range(len(seq) - size + 1)]


def pattern(seq, size):
    return [(i + 1, cg_content(w), information_content(w)) for i, w in sliding_windows(seq, size)]


def center_of_weight(points):
    n = len(points)
    return sum(p[1] for p in points) / n, sum(p[2] for p in points) / n


def read_fasta(path):
    records, name = {}, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                parts = line[1:].split()
                name = parts[3].rstrip(";") if len(parts) > 3 else parts[0]
                records[name] = ""
            elif name:
                records[name] += line.upper()
    return records


def draw_pattern(ax, points, center, title):
    xs = [p[1] for p in points]
    ys = [p[2] for p in points]
    sc = ax.scatter(xs, ys, c=range(len(points)), cmap="viridis", s=14)
    ax.plot(*center, marker="X", markersize=13, color="#c0392b", markeredgecolor="white")
    ax.annotate(f"center ({center[0]:.1f}, {center[1]:.1f})", center, textcoords="offset points",
                xytext=(8, 8), fontsize=8, color="#c0392b")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("C+G %")
    ax.set_ylabel("Information content")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
    return sc


def main():
    p = argparse.ArgumentParser(description="DNA patterns: C+G% and information content per sliding window.")
    p.add_argument("--seq", default=S)
    p.add_argument("--window", type=int, default=30)
    p.add_argument("--promoters", default="promoters.fasta")
    p.add_argument("--promoter", default="GAPDH_1")
    p.add_argument("--no-show", action="store_true")
    a = p.parse_args()
    seq = a.seq.upper()

    print("=" * 70)
    print("DNA PATTERNS: C+G % AND INFORMATION CONTENT")
    print("=" * 70)
    print(f"S = {seq}")
    print(f"Length {len(seq)}, sliding window {a.window}, {len(seq) - a.window + 1} windows")

    cg = cg_content(seq)
    ic = information_content(seq)
    print("\n1. CHECK THE FUNCTIONS ON THE WHOLE S")
    print(f"   C+G %               = {cg:.2f}   (board: {EXPECTED_CG})  " + ("OK" if abs(cg - EXPECTED_CG) < 0.01 else "DIFFERENT"))
    print(f"   Information content = {ic:.2f}   (board: {EXPECTED_IC})  " + ("OK" if abs(ic - EXPECTED_IC) < 0.01 else "DIFFERENT"))
    print(f"   Without the '100 -': {self_match(seq):.2f}")

    pts = pattern(seq, a.window)
    center = center_of_weight(pts)
    print(f"\n2. PATTERN OF S ({len(pts)} windows)")
    print(f"   {'Window':>6}  {'C+G %':>7}  {'IC':>7}   Sequence")
    for pos, c, i in pts:
        print(f"   {pos:>6}  {c:>7.2f}  {i:>7.2f}   {seq[pos - 1:pos - 1 + a.window]}")
    print(f"\n3. CENTER OF WEIGHT OF THE PATTERN: C+G = {center[0]:.2f}, IC = {center[1]:.2f}")

    centers = {"S (lab)": center}
    promoters = read_fasta(a.promoters) if os.path.exists(a.promoters) else {}
    print("\n4. PROMOTER PATTERNS")
    if promoters:
        for name, s in promoters.items():
            c = center_of_weight(pattern(s, a.window))
            centers[name] = c
            print(f"   {name:<9} {len(s)} bp   center C+G {c[0]:6.2f}   IC {c[1]:6.2f}")
    else:
        print(f"   {a.promoters} not found, skipping promoters")

    fig, axes = plt.subplots(1, 2 if a.promoter in promoters else 1, figsize=(12 if a.promoter in promoters else 6, 5.5))
    axes = axes if isinstance(axes, (list, tuple)) or hasattr(axes, "__len__") else [axes]
    draw_pattern(axes[0], pts, center, "Pattern of S (colour = window position)")
    if a.promoter in promoters:
        ps = pattern(promoters[a.promoter], a.window)
        draw_pattern(axes[1], ps, centers[a.promoter], f"Pattern of the {a.promoter.split('_')[0]} promoter")
    plt.tight_layout()
    plt.savefig("pattern.png", dpi=150)

    fig2, ax2 = plt.subplots(figsize=(7, 5.5))
    offsets = {"TP53": (-7, -14, "right"), "HBB": (0, 10, "center"), "BRCA1": (-2, 10, "right"),
               "GAPDH": (2, 10, "left"), "INS": (0, -17, "center"), "ACTB": (9, -4, "left")}
    for name, (x, y) in centers.items():
        label = name.split("_")[0]
        dx, dy, ha = offsets.get(label, (7, 8, "left"))
        ax2.scatter(x, y, s=70, color="#c0392b" if name.startswith("S") else "#1f5fbf")
        ax2.annotate(label, (x, y), textcoords="offset points", xytext=(dx, dy), fontsize=9, ha=ha)
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 100)
    ax2.set_xlabel("C+G % (center)")
    ax2.set_ylabel("Information content (center)")
    ax2.set_title("Centers of weight of all patterns", fontsize=10)
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("centers.png", dpi=150)
    print("\nCharts saved as pattern.png and centers.png")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
