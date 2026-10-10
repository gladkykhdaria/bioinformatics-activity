genetic_code = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",

    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",

    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",

    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

one_letter = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Stop": "*",
}


def to_mrna(seq):
    return seq.upper().replace(" ", "").replace("T", "U")


def read_through(mrna, start):
    letters = ""
    for i in range(start, len(mrna) - 2, 3):
        aa = genetic_code.get(mrna[i:i + 3], "?")
        letters += one_letter.get(aa, "?")
    return letters


def all_proteins(mrna):
    proteins = []
    i = mrna.find("AUG")
    while i != -1:
        protein = []
        j = i
        while j + 3 <= len(mrna):
            aa = genetic_code.get(mrna[j:j + 3], "?")
            if aa == "Stop":
                break
            protein.append(aa)
            j += 3
        proteins.append(protein)
        i = mrna.find("AUG", j + 3)
    return proteins


def first_protein(mrna):
    start = mrna.find("AUG")
    if start == -1:
        return ""
    return read_through(mrna, start).split("*")[0]


def title(text):
    print()
    print("=" * 50)
    print(text)
    print("=" * 50)


def option1_read_through(mrna, start):
    title("OPTION 1: Read-through (mark stops with *, keep going)")
    print("Protein:", read_through(mrna, start))


def option2_fragments(mrna, start):
    title("OPTION 2: Split into fragments at each stop")
    pieces = [p for p in read_through(mrna, start).split("*") if p]
    for n, p in enumerate(pieces, 1):
        print(f"Fragment {n}: {p}")


def option3_all_genes(mrna):
    title("OPTION 3: Find all genes (restart at next AUG after a stop)")
    for n, p in enumerate(all_proteins(mrna), 1):
        letters = "".join(one_letter[a] for a in p)
        print(f"Protein {n}: {'-'.join(p)}  ({letters})")


def option4_count_stops(mrna, start):
    title("OPTION 4: Count stop codons (in frame, from first AUG)")
    counts = {"UAA": 0, "UAG": 0, "UGA": 0}
    for i in range(start, len(mrna) - 2, 3):
        codon = mrna[i:i + 3]
        if codon in counts:
            counts[codon] += 1
    for codon, c in counts.items():
        print(f"{codon}: {c}")
    print("Total:", sum(counts.values()))


def option6_mutation(mrna):
    title("OPTION 6: Detect nonsense mutation (early stop)")
    other = input("Enter mutated sequence (or press Enter to skip): ")
    if not other.strip():
        print("Skipped.")
        return
    normal = first_protein(mrna)
    mutated = first_protein(to_mrna(other))
    print("Normal protein: ", normal, f"({len(normal)} aa)")
    print("Mutated protein:", mutated, f"({len(mutated)} aa)")
    if normal == mutated:
        print("Result: no change in the protein.")
    elif len(mutated) < len(normal) and normal.startswith(mutated):
        print(f"Result: NONSENSE MUTATION - early stop after amino acid {len(mutated)}.")
    else:
        pos = next((k for k in range(min(len(normal), len(mutated)))
                    if normal[k] != mutated[k]), min(len(normal), len(mutated)))
        print(f"Result: protein changed starting at amino acid {pos + 1}.")


def option5_three_frames(mrna):
    title("OPTION 5: Check all 3 reading frames (longest ORF)")
    best = ("", 0)
    for frame in range(3):
        letters = read_through(mrna, frame)
        longest = ""
        for piece in letters.split("*"):
            m = piece.find("M")
            if m != -1 and len(piece) - m > len(longest):
                longest = piece[m:]
        print(f"Frame {frame + 1}: {letters}")
        print(f"   longest ORF: {longest or '-'}")
        if len(longest) > len(best[0]):
            best = (longest, frame + 1)
    if best[0]:
        print(f"Best guess for the gene: frame {best[1]}, {best[0]} ({len(best[0])} aa)")
    else:
        print("No ORF found in any frame.")


mrna = to_mrna(input("Enter DNA or mRNA sequence: "))
start = mrna.find("AUG")

if start == -1:
    print("No start codon (AUG) found.")
else:
    print("mRNA:", mrna)
    option1_read_through(mrna, start)
    option2_fragments(mrna, start)
    option3_all_genes(mrna)
    option4_count_stops(mrna, start)
    option5_three_frames(mrna)
    option6_mutation(mrna)