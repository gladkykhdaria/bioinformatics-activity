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



"""
Gene -> Protein Translator
==========================


  Genome (DNA)       = hard drive    -> permanent storage, stays in the nucleus
  Gene               = a file        -> one stored "program"
  Transcription      = loading to RAM -> DNA is copied into mRNA (T becomes U)
  mRNA               = RAM           -> temporary working copy, used then discarded
  Ribosome           = CPU            -> reads the code and executes it
  Codon (3 bases)    = instruction    -> fixed-size 3-letter word
  tRNA               = data bus       -> fetches the right amino acid for each codon
  AUG                = entry point    -> main(), where execution starts (Met)
  Stop codon         = return / EOF   -> UAA, UAG, UGA end the program
  Protein            = program output -> the chain folds into a 3D working shape

Why 3 letters per codon? 
simply because: 4 bases ^ 2 = 16 combinations, too few for 20 amino
acids. 4 ^ 3 = 64 is enough, so the code has redundancy (6 codons -> Leu),
like several opcodes mapping to the same operation.

How the code mirrors the cell:
  .replace("T", "U")      -> transcription
  find("AUG")             -> ribosome looks for the start signal
  for loop, step 3        -> reading one codon at a time
  genetic_code[codon]     -> tRNA brings the matching amino acid
  protein.append(...)     -> chain grows
  break on "Stop"         -> ribosome releases the finished protein

Life solved information storage and execution about 4 billion years before we did;)

"""


dna = input("Enter DNA or mRNA sequence: ").upper().replace(" ", "")
mrna = dna.replace("T", "U")
start = mrna.find("AUG")

if start == -1:
    print("No start codon (AUG) found.")
else:
    print()
    print(f"{'DNA':<6}{'mRNA':<6}{'Protein':<9}{'1-letter'}")
    print("-" * 29)

    protein = []
    letters = []
    for i in range(start, len(mrna) - 2, 3):
        dna_codon = dna[i:i + 3]
        codon = mrna[i:i + 3]
        amino_acid = genetic_code.get(codon, "?")
        letter = one_letter.get(amino_acid, "?")
        print(f"{dna_codon:<6}{codon:<6}{amino_acid:<9}{letter}")

        if amino_acid == "Stop":
            break
        protein.append(amino_acid)
        letters.append(letter)

    print()
    print("Protein:", "-".join(protein))
    print("1-letter:", "".join(letters))