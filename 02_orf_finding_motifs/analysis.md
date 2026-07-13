# ORF Finding & Motif Scoring

## Overview
Implemented a multi-frame ORF scanner combined with a position weight 
matrix (PWM) motif scorer to evaluate the biological plausibility of 
candidate start codons. ORFs were identified by scanning all 3 reading 
frames for start/stop codon pairs, then scored based on a 13 bp sequence 
context surrounding the start codon.

## Multiple ORFs in a Single Sequence
Several sequences contained multiple candidate ORFs across different, 
overlapping reading frames (e.g. one sequence had 3 ORFs spanning 
positions 46–133, 56–122, and 102–162). Since real genes don't overlap in 
this way across frames, the most likely explanation is chance occurrence 
of start/stop codon pairs — with 3 stop codons out of 64 possible codons 
(~4.7% chance per codon), spurious short ORFs are expected to arise 
randomly, especially in short sequences scanned across all 3 frames.

## Borderline & Failing ORFs
44% of candidate ORFs failed to meet the length and motif score criteria, 
with several scoring in a borderline range (6.5–7.24, just under the 7.25 
cutoff). Sequences scoring just below cutoff (e.g. 7.00) likely represent 
false negatives — a single base substitution in a less critical motif 
position could shift a true ORF just under threshold. ORFs scoring much 
lower are more likely true negatives, suggesting the current threshold may 
be slightly too stringent for borderline cases.

## Statistical Validation

**Probability of a stop codon immediately following a start codon:**  
3/64 ≈ 4.69%

**Probability of no stop codon in 200 codons (600 bp) by chance:**  
(61/64)^200 ≈ 0.0000842 (~0.008%)

**Expected number of stop codons in 200 codons (Poisson λ):**  
λ = 200 × (3/64) = 9.375

**P(0 stop codons) via Poisson:**  
P(X=0) = e^(−9.375) ≈ 0.0000842

The binomial and Poisson approaches converge to the same probability, 
confi