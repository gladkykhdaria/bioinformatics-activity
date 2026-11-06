# Lab 06 - Transposable Elements

## Task

1. Make an artificial DNA sequence of 200 to 400 bases with 3 to 4 transposable elements inside it.
2. Write a program that finds the position (start and end) of each transposable element in that sequence.

## What a transposon is

A transposon (also called a jumping gene) is a piece of DNA that can cut itself out and move to another place in the genome. They were first found by Barbara McClintock in maize.

Two things make a transposon recognisable:

- **Inverted repeats (IR)** at its two ends. The right end is the left end read backwards on the other strand (the reverse complement).
- **Direct repeats** around it. When the transposon is inserted, the host DNA is cut with a few bases between the two strands. After the gaps are filled, these few bases appear on both sides. This is also called the target site duplication (TSD).

```
host DNA:     ...TAG ATGCA CAG...
after:        ...TAG ATGCA [ATGCAGTC .... GACTGCAT] ATGCA CAG...
                     ^^^^^  ^^^^^^^^      ^^^^^^^^  ^^^^^
                     direct inverted      inverted  direct
                     repeat repeat        repeat    repeat
```

Reverse complement of `ATGCAGTC`: reverse it to `CTGACGTA`, then swap A with T and C with G to get `GACTGCAT`.

## Files

| File | What it does |
|---|---|
| `ex1_make_sequence.py` | Builds the sequence. Saves `artificial_dna.fasta` and `answer_key.txt` |
| `ex2_detect_transposons.py` | Reads only the FASTA file and finds the transposons. Uses the answer key only at the end to check |
| `transposon-finder.html` | Both parts in the browser, with a map, coloured sequence and an accuracy test |

## How to run

```
python ex1_make_sequence.py
python ex2_detect_transposons.py
```

Other options:

```
python ex1_make_sequence.py --seed 7 --count 3
python ex2_detect_transposons.py my_sequence.fasta --min-ir 7 --min-tsd 4
```

For the HTML version, open `transposon-finder.html` in a browser. It also has a box where you can paste any sequence.

---

## Exercise 1 - Making the sequence

1. Build 3 or 4 transposons. Each one has a random left IR (8 to 10 bases), a random middle part (15 to 35 bases) and a right IR that is the reverse complement of the left IR.
2. Pick a total length between 200 and 400 and make random host DNA that leaves exactly enough room.
3. Split the host into equal sections and pick one insertion point in each, so the transposons don't overlap.
4. At each point, copy the 5 host bases of the target site, add the transposon, then add the same 5 bases again. This is the direct repeat.
5. Save the sequence and, separately, the true positions.

The answer key is in a separate file on purpose. The detector never reads it while searching, so the check at the end is fair.

## Exercise 2 - Detecting them

1. **Index.** Save every 7-letter word of the sequence in a dictionary, with its positions.
2. **Find IR pairs.** For each word, compute its reverse complement and look it up. A match 24 to 150 bases to the right is a possible transposon.
3. **Grow.** Extend the match outwards and inwards while the bases are still reverse complements. This gives the full IR.
4. **Check the direct repeats.** The 4 to 10 bases right before the start must be the same as the bases right after the end. If they are not, move both edges one base inwards and try again.
5. **Remove overlaps.** If two candidates overlap, keep the one with the higher score:

$$
\text{score} = \text{IR length} + 2 \times \text{direct repeat length}
$$

6. Print start and end, counting from 1.

Step 4 is needed because sometimes the base just outside the IR happens to match too, and the IR looks one base longer than it really is. At that wrong edge the direct repeats don't match, so the program moves in by one and finds the right edge.

## Results

Default run (seed 6, 4 transposons, 341 bases):

| # | Start | End | Length | Left IR | Right IR | Direct repeat |
|---|---|---|---|---|---|---|
| 1 | 30 | 83 | 54 | AGCTGGAACA | TGTTCCAGCT | GGGTT |
| 2 | 114 | 163 | 50 | TACGTTCAAA | TTTGAACGTA | TTACA |
| 3 | 211 | 253 | 43 | ATGAACTGGA | TCCAGTTCAT | AAAGA |
| 4 | 307 | 341 | 35 | CAAGGTAC | GTACCTTG | ATATG |

All 4 match the answer key exactly, with 0 false hits.

Small detail: for #2 the detector reports a 10-base IR, but the real one was 9. The next base in the middle part happened to fit the pattern too. The start and end are still exactly right, because only the outer edge matters for the position.

### Testing it on many sequences

I generated 14,000 transposons in 4,000 random sequences and ran the detector on all of them:

| | |
|---|---|
| Found with exact start and end | 13,990 of 14,000 (99.93%) |
| False hits | 10 |

In the few misses, the random bases right next to a transposon happened to look like both an extra piece of IR and a direct repeat at the same time, so the edges came out 2 bases too wide.

## Why random DNA almost never looks like a transposon

A 7-base word can be written in $4^7 = 16{,}384$ ways. In a 400-base sequence, looking about 125 bases ahead, the number of random IR pairs we expect is:

$$
\frac{400 \times 125}{16{,}384} \approx 3
$$

Each of these also needs a direct repeat of at least 4 bases in exactly the right place. The chance of that is $\frac{1}{4^4} = \frac{1}{256}$. So the expected number of false transposons per sequence is:

$$
\frac{3}{256} \approx 0.01
$$

That fits the test: 10 false hits in 4,000 sequences is 0.0025 per sequence, which is even lower.

## What I noticed

- One signal alone is not enough. Random 7-base inverted repeats appear a few times in every sequence. Only the combination with direct repeats makes the detection reliable.
- The direct repeat is the most useful clue for finding the exact edges.
- The dictionary makes the search fast. Detection on one sequence takes under 1 ms.

## Limitations

- My IRs are perfect copies. In real genomes old transposons collect mutations, so the IRs are not exact anymore.
- I only simulate DNA transposons (cut and paste). Retrotransposons, which copy themselves through RNA, look different and don't have inverted repeats at the ends.
- My minimum direct repeat is 4 bases. Some real transposon families make only a 2-base duplication, and my detector would miss them.

## Questions I had for the professor

1. In real genomes the inverted repeats of old transposons have mutations. How many mismatches should a detector allow before random DNA starts to look like transposons too?
2. Most transposable elements in the human genome are retrotransposons, which have long direct repeats (LTRs) at their ends or no repeats at all. Does that mean detection based on inverted repeats only finds a small part of them?
3. I read that different transposon families make direct repeats of a typical length, for example only `TA` (2 bases) for the mariner family. Can the length of the direct repeat be used to tell which family a transposon belongs to?
4. If I lower my minimum direct repeat to 2 bases to catch families like mariner, false hits go up a lot. How do real tools handle this trade-off?
