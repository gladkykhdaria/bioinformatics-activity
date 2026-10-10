import random

alphabet = "ACGT"
S = "".join(random.choice(alphabet) for _ in range(20))
print("S =", S)

counts = {letter: 0 for letter in alphabet}
for ch in S:
    counts[ch] += 1

for letter in alphabet:
    f = counts[letter] / len(S)
    print(f"{letter}: count={counts[letter]}, freq={f:.2f} ({f*100:.0f}%)")