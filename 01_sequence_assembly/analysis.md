# Sequence Assembly: Greedy Overlap-Layout Algorithm

## Overview
I built a greedy overlap-layout assembler from scratch in Python. The 
algorithm works by starting with an unused read, then repeatedly checking all remaining reads to find the best overlap (either extending left or right), and merging that read in. This continues until every read has been incorporated into a contig. I also calculated read depth and density across the assembled contigs to visualize coverage.

## Read Density Analysis
When I plotted read density across the contig, I saw lower coverage at the 
ends (around positions 1-20 and 100-130) and higher, more consistent 
coverage in the middle (around positions 20-100). This isn't sequencing 
bias, it's expected in assembled contigs. Reads that cover the ends of a 
contig are constrained to start or end right at those positions, but reads 
covering the middle can start anywhere before or after that point, so more 
reads end up overlapping the middle by chance.

## K-Parameter Sensitivity
I tested how changing k (the minimum overlap length) affected the number 
of contigs produced:

| k value | Number of contigs |
|---|---|
| 10 | 20 |
| 9 | 20 |
| 8 | 20 |
| 7 | 20 |
| 6 | 20 |
| 5 | 20 |
| 4 | 19 |
| 3 | 13 |

| k value | Number of contigs |
|---|---|
| 10 | 20 |
| 11 | 20 |
| 12 | 45 |
| 13 | 63 |
| 14 | 63 |

The contig count stayed the same from k=7 through k=11, but changed a lot 
at both extremes:

- **Lowering k** (k=3-4): fewer contigs, because short overlaps start 
  matching by chance rather than because reads actually came from the same 
  region. At k=3 especially, common motifs like the start codon ATG can 
  create false overlaps just by random chance.
- **Raising k** (k=12-14): way more contigs, because requiring a longer 
  overlap means some real overlaps get missed and reads that should've 
  merged stay separate.

This shows the trade-off with k: it needs to be small enough to catch real 
overlaps, but large enough to avoid matching sequences that just happen to 
share a few bases. How high you can push k really depends on sequencing 
depth and coverage since with more coverage, you're more likely to get long 
true overlaps, so you can use a higher k without missing anything.

## Paired-End vs. Mate-Pair Sequencing

Given a 400 bp fragment, 80 bp reads, and a 3,000 bp mate-pair library:

- **Paired-end**: L = 400 − 80 − 80 = **240 bp**
- **Mate-pair**: L = 3,000 − 80 − 80 = **2,840 bp**

Mate-pair sequencing gets a much bigger insert size because it works 
differently by circularizing a large DNA fragment first, then sequencing 
from the joined ends, so it's not limited by fragment size the way regular 
paired-end prep is. That larger distance is really useful for resolving 
long repetitive regions during assembly, since you can span across a 
repeat that a shorter paired-end read wouldn't be able to.

## Note
Dataset is simulated for coursework purposes.