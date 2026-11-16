import argparse
import gzip
import math

import matplotlib.pyplot as plt

DEFAULT_GENOME = "../genomes/NC_005816.fasta"

ENZYMES = [
    {"name": "EcoRI", "site": "GAATTC", "top": 1, "bottom": 5, "source": "Escherichia coli RY13"},
    {"name": "BamHI", "site": "GGATCC", "top": 1, "bottom": 5, "source": "Bacillus amyloliquefaciens H"},
    {"name": "HindIII", "site": "AAGCTT", "top": 1, "bottom": 5, "source": "Haemophilus influenzae Rd"},
    {"name": "TaqI", "site": "TCGA", "top": 1, "bottom": 3, "source": "Thermus aquaticus"},
    {"name": "HaeIII", "site": "GGCC", "top": 2, "bottom": 2, "source": "Haemophilus aegyptius"},
]

LADDER = [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000]
GEL_MIN, GEL_MAX = 100, 12000

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def revcomp(s):
    return "".join(COMP[c] for c in reversed(s))


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
    seq = "".join(parts)
    return name, "".join(c for c in seq if c in "ACGT")


def find_sites(dna, site, circular):
    text = dna + dna[:len(site) - 1] if circular else dna
    hits, i = [], text.find(site)
    while i != -1 and i < len(dna):
        hits.append(i)
        i = text.find(site, i + 1)
    return hits


def cut_positions(sites, offset, length):
    return sorted({(s + offset) % length for s in sites})


def fragments(cuts, length, circular):
    if not cuts:
        return [length]
    if circular:
        return [(cuts[(k + 1) % len(cuts)] - cuts[k]) % length or length for k in range(len(cuts))]
    edges = [0] + cuts + [length]
    return [b - a for a, b in zip(edges, edges[1:]) if b > a]


def expected_sites(dna, site):
    p = {b: dna.count(b) / len(dna) for b in "ACGT"}
    return len(dna) * math.prod(p[b] for b in site)


def end_type(e):
    if e["top"] == e["bottom"]:
        return "blunt", ""
    over = e["site"][min(e["top"], e["bottom"]):max(e["top"], e["bottom"])]
    return ("5' sticky" if e["top"] < e["bottom"] else "3' sticky"), over


def digest(dna, enzymes, circular):
    owner = {}
    for e in enzymes:
        for c in cut_positions(find_sites(dna, e["site"], circular), e["top"], len(dna)):
            owner.setdefault(c, []).append(e["name"])
    cuts = sorted(owner)
    return cuts, fragments(cuts, len(dna), circular), owner


def print_mix(dna, enzymes, circular):
    cuts, frags, owner = digest(dna, enzymes, circular)
    label = " + ".join(e["name"] for e in enzymes)
    print(f"\n{label}")
    print(f"  {len(cuts)} cuts, {len(frags)} fragments")
    print("  cuts (enzyme): " + ", ".join(f"{c:,} ({'/'.join(owner[c])})" for c in cuts))
    print("  fragments (bp): " + ", ".join(f"{f:,}" for f in sorted(frags, reverse=True)))
    print(f"  check: sum = {sum(frags):,} bp")
    return label, frags


def migration(length):
    lo, hi = math.log10(GEL_MIN), math.log10(GEL_MAX)
    return 4 + 92 * (hi - math.log10(length)) / (hi - lo)


def draw_gel(lanes, path):
    n = len(lanes) + 1
    fig, ax = plt.subplots(figsize=(1.05 * n + 1.6, 7.5))
    ax.set_facecolor("#0f1d2e")
    for k in range(n):
        ax.add_patch(plt.Rectangle((k + 0.12, -3.2), 0.76, 2.2, color="#24364d"))
    for size in LADDER:
        y = migration(size)
        ax.add_patch(plt.Rectangle((0.12, y - 0.35), 0.76, 0.7, color="#9fd3ff", alpha=0.85))
        ax.text(-0.12, y, f"{size:,}", color="#b8c4d4", ha="right", va="center", fontsize=7.5)
    for k, (label, frags) in enumerate(lanes, start=1):
        counts = {}
        for f in frags:
            counts[f] = counts.get(f, 0) + 1
        for size, c in counts.items():
            if size < GEL_MIN:
                continue
            y = migration(min(size, GEL_MAX * 0.98))
            alpha = min(1.0, 0.25 + 0.75 * (size * c) / 4000)
            ax.add_patch(plt.Rectangle((k + 0.12, y - 0.45), 0.76, 0.9, color="#7dff9e", alpha=alpha))
        small = sum(1 for f in frags if f < GEL_MIN)
        if small:
            ax.text(k + 0.5, 99.5, f"+{small} ran off", color="#9fb3c8", ha="center", fontsize=6.5)
    ax.set_xlim(-1.0, n)
    ax.set_ylim(101, -4.5)
    ax.set_xticks([k + 0.5 for k in range(n)])
    ax.set_xticklabels(["Ladder"] + [l for l, _ in lanes], fontsize=8.5, rotation=30, ha="right")
    ax.set_yticks([])
    ax.text(-0.55, -2.1, "bp", color="#b8c4d4", ha="center", fontsize=8)
    ax.set_title("Simulated agarose gel of the restriction digests", fontsize=11)
    for s in ax.spines.values():
        s.set_visible(False)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    return fig


def main():
    p = argparse.ArgumentParser(description="Digest a DNA sequence with 5 restriction enzymes and simulate the gel.")
    p.add_argument("genome", nargs="?", default=DEFAULT_GENOME)
    p.add_argument("--linear", action="store_true")
    p.add_argument("--mix", action="append", default=[], metavar="ENZ1,ENZ2")
    p.add_argument("--no-show", action="store_true")
    a = p.parse_args()

    name, dna = read_fasta(a.genome)
    circular = not a.linear
    L = len(dna)

    print("=" * 72)
    print("RESTRICTION DIGEST")
    print("=" * 72)
    print("DNA:     ", name)
    print("Length:  ", f"{L:,} bp", "(circular)" if circular else "(linear)")

    print("\nENZYMES")
    print(f"{'Name':<9}{'Site':<8}{'Cut':<12}{'Ends':<11}{'Overhang':<9}{'Palindrome':<11}Comes from")
    print("-" * 72)
    for e in ENZYMES:
        kind, over = end_type(e)
        shown = e["site"][:e["top"]] + "^" + e["site"][e["top"]:]
        pal = "yes" if revcomp(e["site"]) == e["site"] else "NO"
        print(f"{e['name']:<9}{e['site']:<8}{shown:<12}{kind:<11}{over or '-':<9}{pal:<11}{e['source']}")

    print("\n1. NUMBER OF CLEAVAGES")
    print(f"{'Enzyme':<9}{'Cuts':>5}{'Expected':>10}   Fragments")
    print("-" * 72)
    lanes = []
    for e in ENZYMES:
        sites = find_sites(dna, e["site"], circular)
        cuts = cut_positions(sites, e["top"], L)
        frags = fragments(cuts, L, circular)
        e["sites"], e["cuts"], e["frags"] = sites, cuts, frags
        lanes.append((e["name"], frags))
        print(f"{e['name']:<9}{len(cuts):>5}{expected_sites(dna, e['site']):>10.1f}   {len(frags)}")

    print("\n2. CLEAVAGE POSITIONS AND FRAGMENT LENGTHS")
    for e in ENZYMES:
        print(f"\n{e['name']}  ({e['site']}, cuts after base {e['top']} of the site)")
        if not e["cuts"]:
            print("  no sites: the DNA stays in one piece")
            continue
        print("  cut after base: " + ", ".join(f"{c:,}" for c in e["cuts"]))
        print("  fragments (bp): " + ", ".join(f"{f:,}" for f in sorted(e["frags"], reverse=True)))
        print(f"  check: sum = {sum(e['frags']):,} bp")

    print("\n3. SEVERAL ENZYMES AT ONCE")
    by_name = {e["name"].lower(): e for e in ENZYMES}
    mixes = a.mix or ["EcoRI,BamHI", "EcoRI,HindIII", "BamHI,HindIII"]
    for m in mixes:
        chosen = [by_name[n.strip().lower()] for n in m.split(",") if n.strip().lower() in by_name]
        if len(chosen) >= 2:
            lanes.append(print_mix(dna, chosen, circular))
        else:
            print(f"\n  skipped '{m}': give at least 2 of {', '.join(e['name'] for e in ENZYMES)}")
    all_cuts, all_frags, _ = digest(dna, ENZYMES, circular)
    lanes.append(("All 5", all_frags))
    print(f"\nAll 5 enzymes together: {len(all_cuts)} cuts, {len(all_frags)} fragments, "
          f"largest {max(all_frags):,} bp, smallest {min(all_frags)} bp")

    print("\n4. GEL")
    for label, frags in lanes:
        seen = sorted(set(f for f in frags if f >= GEL_MIN), reverse=True)
        lost = sum(1 for f in frags if f < GEL_MIN)
        print(f"  {label:<8} {len(seen):>2} visible bands" + (f", {lost} too small (< {GEL_MIN} bp)" if lost else ""))
    draw_gel(lanes, "gel_digest.png")
    print("\nGel saved as gel_digest.png")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
