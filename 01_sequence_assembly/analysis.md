# Sequence Assembly: Greedy Overlap-Layout Algorithm

## Overview
Implemented a greedy overlap-layout assembler from scratch in Python. The 
algorithm iteratively extends contigs by finding the best-scoring overlap 
(left or right) among remaining unused reads, continuing until no reads 
remain. Read depth and density were also calculated across assembled 
contigs to visualize coverage patterns.

## Read Density Analysis
Read density across the contig showed lower coverage at the ends 
(positions ~1-20 and ~100-130) and higher, more consistent coverage in the 
middle (positions ~20-100). This pattern is expected in assembled contigs 
rather than indicative of sequencing bias — reads covering the ends of a 
contig are constrained to start or end at those positions, while reads 
covering the middle can originate from a wider range of start positions, 
naturally inflating middle coverage.

## K-Parameter Sensitivity

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

Contig count remained stable across a middle range of k (7–11), but changed 
sharply at the extremes:

- **Decreasing k** (k=3-4): Contig count dropped as short, low-specificity 
  overlaps caused reads to merge that likely didn't originate from the same 
  genomic region. At k=3, common short motifs (e.g. the start codon ATG) 
  create frequent false-positive overlaps by chance.
- **Increasing k** (k=12-14): Contig count rose sharply as the overlap 
  requirement became too strict to detect true overlaps, fragmenting the 
  assembly.

This reflects a fundamental trade-off in overlap-based assembly: k must be 
low enough to detect real overlaps but high enough to avoid spurious 
matches. The appropriate upper bound for k is directly tied to sequencing 
depth and coverage — higher coverage supports higher k, since true overlaps 
are more likely to be long, while lower coverage requires a lower k to 
avoid missing real overlaps.

## Paired-End vs. Mate-Pair Sequencing

Given a 400 bp fragment, 80 bp reads, and a 3,000 bp mate-pair library:

- **Paired-end**: L = 400 − 80 − 80 = **240 bp**
- **Mate-pair**: L = 3,000 − 80 − 80 = **2,840 bp**

Mate-pair sequencing achieves much larger insert sizes because it relies on 
circularizing large DNA fragments before sequencing the joined ends, rather 
than being limited by fragment size as in standard paired-end prep. This 
allows mate-pair reads to span much larger genomic distances, which is 
useful for resolving long repetitive regions during assembly.

## Note
Dataset is simulated for coursework purposes.