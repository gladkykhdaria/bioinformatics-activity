import argparse

import matplotlib.pyplot as plt

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"


def score(a, b, match, mismatch):
    return match if a == b else mismatch


def init_matrix(n, m, gap):
    M = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        M[i][0] = i * gap
    for j in range(m + 1):
        M[0][j] = j * gap
    return M


def fill(s1, s2, match, mismatch, gap):
    n, m = len(s1), len(s2)
    M = init_matrix(n, m, gap)
    ways = [[1 if i == 0 or j == 0 else 0 for j in range(m + 1)] for i in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = M[i - 1][j - 1] + score(s1[i - 1], s2[j - 1], match, mismatch)
            up = M[i - 1][j] + gap
            left = M[i][j - 1] + gap
            best = max(diag, up, left)
            M[i][j] = best
            ways[i][j] = ((ways[i - 1][j - 1] if diag == best else 0)
                          + (ways[i - 1][j] if up == best else 0)
                          + (ways[i][j - 1] if left == best else 0))
    return M, ways


def traceback(s1, s2, M, match, mismatch, gap):
    i, j = len(s1), len(s2)
    top, bottom, path = [], [], [(i, j)]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and M[i][j] == M[i - 1][j - 1] + score(s1[i - 1], s2[j - 1], match, mismatch):
            top.append(s1[i - 1])
            bottom.append(s2[j - 1])
            i, j = i - 1, j - 1
        elif i > 0 and M[i][j] == M[i - 1][j] + gap:
            top.append(s1[i - 1])
            bottom.append("-")
            i -= 1
        else:
            top.append("-")
            bottom.append(s2[j - 1])
            j -= 1
        path.append((i, j))
    return "".join(reversed(top)), "".join(reversed(bottom)), list(reversed(path))


def middle_line(a, b):
    return "".join("|" if x == y else (" " if "-" in (x, y) else ".") for x, y in zip(a, b))


def print_matrix(s1, s2, M, path):
    on_path = set(path)
    print("       " + "     -" + "".join(f"{c:>6}" for c in s2))
    for i in range(len(s1) + 1):
        label = "-" if i == 0 else s1[i - 1]
        row = "".join((f"[{M[i][j]:>3}]" if (i, j) in on_path else f" {M[i][j]:>3} ").rjust(6) for j in range(len(s2) + 1))
        print(f"   {label:>2} " + row)


def draw(s1, s2, M, path, a, b, file):
    n, m = len(s1), len(s2)
    fig, ax = plt.subplots(figsize=(0.48 * (m + 1) + 2, 0.48 * (n + 1) + 2))
    ax.imshow(M, cmap="Blues", aspect="equal")
    for i in range(n + 1):
        for j in range(m + 1):
            ax.text(j, i, M[i][j], ha="center", va="center", fontsize=7.5,
                    color="white" if M[i][j] > max(map(max, M)) * 0.6 else "#1b2a3a")
    ys, xs = zip(*path)
    ax.plot(xs, ys, color="#1e8a4c", linewidth=2.5, marker="o", markersize=4)
    ax.set_xticks(range(m + 1))
    ax.set_xticklabels(["-"] + list(s2), fontsize=9)
    ax.set_yticks(range(n + 1))
    ax.set_yticklabels(["-"] + list(s1), fontsize=9)
    ax.xaxis.tick_top()
    ax.set_xlabel("S2")
    ax.set_ylabel("S1")
    ax.set_title(f"Needleman-Wunsch score matrix and traceback\n{a}\n{middle_line(a, b)}\n{b}",
                 fontsize=9, family="monospace", pad=28)
    plt.tight_layout()
    plt.savefig(file, dpi=150)
    return fig


def main():
    p = argparse.ArgumentParser(description="Global alignment of two DNA sequences with Needleman-Wunsch.")
    p.add_argument("--s1", default=S1)
    p.add_argument("--s2", default=S2)
    p.add_argument("--match", type=int, default=1)
    p.add_argument("--mismatch", type=int, default=-1)
    p.add_argument("--gap", type=int, default=-1)
    p.add_argument("--no-show", action="store_true")
    a = p.parse_args()

    s1, s2 = a.s1.upper(), a.s2.upper()
    M, ways = fill(s1, s2, a.match, a.mismatch, a.gap)
    top, bottom, path = traceback(s1, s2, M, a.match, a.mismatch, a.gap)
    mid = middle_line(top, bottom)

    print("=" * 70)
    print("NEEDLEMAN-WUNSCH GLOBAL ALIGNMENT")
    print("=" * 70)
    print(f"S1 = {s1}  ({len(s1)} bases)")
    print(f"S2 = {s2}  ({len(s2)} bases)")
    print(f"Scores: match {a.match:+d}, mismatch {a.mismatch:+d}, gap {a.gap:+d}")

    print("\n1. SCORE MATRIX  ([x] = cell on the traceback path)\n")
    print_matrix(s1, s2, M, path)

    print("\n2. TRACEBACK PATH (row, column), from the top-left corner to the bottom-right")
    print("   " + " -> ".join(f"({i},{j})" for i, j in path))

    matches = mid.count("|")
    mismatches = mid.count(".")
    gaps = mid.count(" ")
    print("\n3. ALIGNMENT\n")
    print("   S1  " + top)
    print("       " + mid)
    print("   S2  " + bottom)
    print(f"\n   Score:       {M[-1][-1]}")
    print(f"   Matches:     {matches}")
    print(f"   Mismatches:  {mismatches}")
    print(f"   Gaps:        {gaps}")
    print(f"   Identity:    {matches}/{len(top)} = {matches / len(top) * 100:.1f}%")
    print(f"   Equally good alignments with this score: {ways[-1][-1]}")

    draw(s1, s2, M, path, top, bottom, "alignment_matrix.png")
    print("\nMatrix picture saved as alignment_matrix.png")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
