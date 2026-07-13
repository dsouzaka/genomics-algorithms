# ORF Finding & Motif Scoring

## Overview
I built a multi-frame ORF scanner combined with a position weight matrix 
(PWM) motif scorer to check how biologically plausible each candidate 
start codon actually is. The scanner looks across all 3 reading frames for 
start/stop codon pairs, then scores each one based on the 13 bp sequence 
context around the start codon.

## Multiple ORFs in a Single Sequence
A few sequences had multiple candidate ORFs across different, overlapping 
reading frames and one sequence even had 3 ORFs spanning positions 46–133, 
56–122, and 102–162. Since real genes don't overlap like that across 
frames, I think the most likely explanation is that these are just chance 
occurrences of start/stop codon pairs. With 3 stop codons out of 64 
possible codons (about a 4.7% chance per codon), you'd expect some short, 
spurious ORFs to show up randomly, especially since I'm scanning short 
sequences across all 3 frames.

## Borderline & Failing ORFs
44% of candidate ORFs didn't meet the length and motif score criteria, and 
a handful scored in a borderline range (6.5–7.24, just under the 7.25 
cutoff). I think the ones scoring just below cutoff (like 7.00) are likely 
false negatives. A single base substitution in a less critical motif 
position could easily push a real ORF just under the threshold. The ones 
scoring much lower are probably true negatives, which makes me think the 
current threshold might be a little too strict for the borderline cases.

## Statistical Validation

**Probability of a stop codon immediately following a start codon:**  
3/64 ≈ 4.69%

**Probability of no stop codon in 200 codons (600 bp) by chance:**  
(61/64)^200 ≈ 0.0000842 (~0.008%)

**Expected number of stop codons in 200 codons (Poisson λ):**  
λ = 200 × (3/64) = 9.375

**P(0 stop codons) via Poisson:**  
P(X=0) = e^(−9.375) ≈ 0.0000842

The binomial and Poisson approaches give basically the same probability, 
which confirms that a long ORF (600 bp) showing up by chance in a random 
sequence would be extremely unlikely. This shows the filtering criteria are 
actually doing a good job separating real ORFs from random noise.

## Simulated ORF Length Distribution
I generated 1,000 random 150 bp sequences and looked at the ORF length 
distribution (1,180 ORFs total, mean length 38.4 bp). The distribution was 
strongly right-skewed, which matches the ~4.69% per-codon stop 
probability since most random ORFs terminate quickly, and longer ones get 
rarer and rarer.

Only 22.5% of the random 150 bp sequences had an ORF longer than 60 bp, 
compared to 70% of the real contigs passing both the length and motif 
score thresholds in the actual dataset. That gap (22.5% expected by chance 
vs. 70% observed) is a good sign that the ORFs I found in the real data 
are genuine biological signal, not just random noise.

## Note
Dataset is simulated for coursework purposes.