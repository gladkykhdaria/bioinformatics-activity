import random
from itertools import product
import matplotlib.pyplot as plt

alphabet = "ACGT"
S = ""
for _ in range(20):
    S += random.choice(alphabet)

print("S =", S)

def kmer_freq(seq, k):
    combos = ["".join(p) for p in product(alphabet, repeat=k)]
    counts = {c: 0 for c in combos}
    for i in range(len(seq) - k + 1):
        counts[seq[i:i + k]] += 1
    total = len(seq) - k + 1
    return {c: counts[c] / total for c in combos}, counts

freq2, cnt2 = kmer_freq(S, 2)
freq3, cnt3 = kmer_freq(S, 3)

print("\n2-letter combinations (16):")
for c in freq2:
    print(f"{c}: count={cnt2[c]}, freq={freq2[c]:.3f}")

print("\n3-letter combinations (64):")
for c in freq3:
    print(f"{c}: count={cnt3[c]}, freq={freq3[c]:.3f}")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))
ax2.tick_params(axis="x", rotation=90, labelsize=6)

ax1.bar(freq2.keys(), freq2.values(), color="steelblue")
ax1.set_title(f"Dinucleotide frequencies (S = {S})")
ax1.set_ylabel("Frequency")

ax2.bar(freq3.keys(), freq3.values(), color="seagreen")
ax2.set_title("Trinucleotide frequencies")
ax2.set_ylabel("Frequency")
ax2.tick_params(axis="x", rotation=90, labelsize=8)

plt.tight_layout()
plt.show()