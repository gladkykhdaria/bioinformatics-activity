# Lab 05 - Gel Electrophoresis and DNA Repeats

## Sequence used

*Ageratum yellow vein China virus*, NCBI RefSeq NC_021202.1, 2,739 nt (the task asks for 1000 to 3000 nt).

## Files

| File | What it does |
|---|---|
| `ex1_gel_electrophoresis.py` | Exercise 1: takes 10 random fragments and draws the gel (`gel.png`) |
| `ex2_tandem_repeats.py` | Exercise 2: finds repeats with a 3 to 5 base unit |
| `gel-repeats-lab.html` | Both exercises in the browser, with an animated gel and the repeats marked in the sequence |
| `gel.png` | Gel picture made by the Python script |

## How to run

```
python ex1_gel_electrophoresis.py
python ex2_tandem_repeats.py
```

Both use the lab sequence from `../genomes/` by default. You can also give another file:

```
python ex1_gel_electrophoresis.py ../genomes/other.fna.gz --seed 7
python ex2_tandem_repeats.py ../genomes/other.fna.gz --min 6 --max 10
```

For the HTML version, open `gel-repeats-lab.html` in a browser.

---

## Exercise 1 - Gel electrophoresis

### How I did it

1. Read the sequence from the FASTA file.
2. For each of the 10 samples, pick a random length between 100 and 3000 (but not more than 2,739, the length of the sequence) and a random start position. Copy that part.
3. Store all 10 fragments in a list.
4. Calculate how far each fragment travels in the gel.
5. Draw the gel with a ladder lane next to the samples.

### How the movement works

DNA has a negative charge, so in an electric field it moves towards the + side. The gel is like a net. Short pieces get through it easily and travel far. Long pieces get stuck more often and stay closer to the well.

The distance does not change in a straight line with the size. It follows the logarithm of the size:

$$
d = a - b \cdot \log_{10}(L)
$$

where $L$ is the fragment length in bp. I chose $a$ and $b$ so that 100 bp ends near the bottom of the gel and 3000 bp stays near the top.

This means 100 bp and 200 bp are far apart on the gel, but 2000 bp and 2100 bp are almost at the same place.

### Result (Python, seed 42)

The array: `[2719, 202, 1103, 671, 2333, 2518, 230, 483, 1052, 2565]`

Order on the gel from top to bottom:
2719 > 2565 > 2518 > 2333 > 1103 > 1052 > 671 > 483 > 230 > 202

The four longest fragments (2333 to 2719) are very close to each other at the top. The two shortest (202 and 230) are low and easy to see. This matches what happens on a real gel.

(The HTML version uses a different random generator, so with the same seed it gives different fragments. The logic is the same.)

---

## Exercise 2 - Repeats

### What I looked for

A short unit of 3, 4 or 5 bases written two or more times in a row, so that the whole repeated block is 6 to 10 bases long. Examples:

| Example | Unit | Copies | Length |
|---|---|---|---|
| ATG ATG | ATG | 2 | 6 |
| ATG ATG ATG | ATG | 3 | 9 |
| GATC GATC | GATC | 2 | 8 |
| ACGTT ACGTT | ACGTT | 2 | 10 |

### How I did it

1. For each unit length $k$ (3, 4 and 5), go through the sequence and compare every letter with the letter $k$ positions later.
2. If they are the same many times in a row, the sequence is repeating with that unit. The repeated region is the number of matches plus $k$, and the number of full copies is the region length divided by $k$.
3. Keep it only if there are at least 2 copies and the length is between 6 and 10.
4. Skip units that are made of an even smaller repeat, like `AAA` (just A) or `ATAT` (just AT).
5. Compare the result with random sequences that have the same letters.

I did it this way and not by searching every 3-letter word twice in a row, because that would count the same region more than once. For example `ATGATGA` contains both `ATG ATG` and `TGA TGA`.

### Result

| | |
|---|---|
| Repeats found | 46 |
| 3-base units | 37 |
| 4-base units | 9 |
| 5-base units | 0 |
| Bases inside repeats | 291 of 2,739 (10.6%) |
| Most common unit | GAA (5 times) |
| Repeats longer than 10 bases | none |

Some examples: `CAC CAC CAC` at position 32, `ATG ATG ATG` at 1672, `GAA GAA GAA` at 2040, `ATAC ATAC` at 220.

### Compared with random DNA

I shuffled the same letters into a random order 100 times and counted the repeats again. The average was about 42, with a range of 28 to 55. The virus has 46, which is inside this range. So most of these short repeats would also appear just by luck. This virus does not have special repeat regions.

It also makes sense that 3-base units are the most common. The chance that 3 letters repeat right after themselves is $(1/4)^3$. For 5 letters it is $(1/4)^5$, which is 16 times smaller.

---

## What I noticed

- In Exercise 1, fragments with very different sizes like 2333 and 2719 end up very close on the gel. A long gel or a lower agarose concentration would be needed to tell them apart well.
- In the HTML version, two of my fragments were 2350 and 2384 bp, only 34 bp apart. On my simulated gel they almost touch. On a real gel they would probably look like one band.
- In Exercise 2, the number of repeats in the virus is almost the same as in random DNA, so the repeats are not doing anything special here.

## Questions I had for the professor

1. My model says distance depends on $\log_{10}$ of the length. I read that for very long DNA this stops working and all fragments end up together near the top. Is this why pulsed-field gel electrophoresis is used for very large DNA?
2. My fragments are cut at random places. In a real lab they would come from restriction enzymes that cut at a specific sequence, like `GAATTC` for EcoRI. Would it make sense to find these sites in the sequence and predict the gel from them?
3. In my results, 3-base repeats are much more common than 4 or 5. Is this only because they are more likely by chance, or also because inside genes a 3-base repeat does not shift the reading frame, so it is less harmful?
4. Some diseases, like Huntington's disease, are caused by a 3-base repeat (`CAG`) that becomes too long. How many copies are normal, and at what point does it become a problem?
