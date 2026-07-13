# Phylogenetics & Molecular Evolution

## Overview
This module covers two related techniques: UPGMA clustering to build a 
phylogenetic tree from pairwise distances, and a codon-based log-odds 
scoring model to classify sequences as coding or non-coding based on 
ancestor-to-descendant mutation patterns.

## UPGMA Phylogenetic Tree

Given a 6-species distance matrix, UPGMA clustering produced the following 
tree (Newick format):
(M_Spacii:6.875,(((T_Pain:1.0,G_Unit:1.0):2.0,Q_Doba:3.0):1.125,(R_Mani:2.0,A_Finch:2.0):2.125):2.75)
![UPGMA Tree](upgma_tree.png)

The algorithm iteratively merges the two closest clusters (by average 
pairwise distance) until a single root remains, producing an ultrametric 
tree where all leaves are equidistant from the root. T_Pain and G_Unit 
clustered first as the most closely related pair, followed by progressive 
merging with more distantly related species/clusters.

## Coding vs. Non-Coding Classification

Implemented log-odds scoring to classify ancestor-to-M.-Spacii codon 
transitions as coding or non-coding, using pre-trained 64x64 codon 
transition matrices for each model.

### Misclassification with Increased Divergence

Testing on the original dataset (`spacii.fa`) produced accurate 
classifications, but testing on a more divergent dataset 
(`spacii_2100.fa`) resulted in 100% of sequences — including true coding 
sequences — being classified as non-coding.

This happens because:
- Coding scores for the original dataset ranged from about -50 to -85, 
  while the divergent dataset's scores ranged from about -300 to -430
- More negative scores indicate the coding model considers the observed 
  mutations increasingly unlikely to have occurred in a coding region
- Coding regions are under purifying selection (mutations affecting amino 
  acid structure are selected against), but given enough evolutionary time, 
  even coding sequences accumulate mutations that begin to resemble 
  random, non-coding-like patterns
- The model was calibrated for a specific level of divergence; applying it 
  to more divergent sequences means the probability estimates no longer 
  reflect the actual evolutionary distance, causing systematic 
  misclassification

### ROC Curve & Threshold Selection

![ROC Curve](roc_curve.png)

The ROC curve shows the trade-off between true positive rate (TPR) and 
false positive rate (FPR) as the classification threshold changes:

| Threshold | TPR | FPR |
|---|---|---|
| -99.69 | 0.684 | 0.048 |
| -108.9 | 1.0 | 0.19 |

The threshold at -99.69 represents a good balance — the "elbow" of the 
curve, where TPR is maximized before FPR begins increasing sharply. At 
this point, about 68% of coding sequences are correctly identified with 
only ~5% of non-coding sequences misclassified.

Lowering the threshold to -108.9 achieves 100% sensitivity (no true coding 
sequences missed) at the cost of a higher 19% false positive rate. Given 
that missing a real coding gene (false negative) means losing a functional 
component from analysis entirely, while a false positive only costs extra 
downstream analysis time, **sensitivity was prioritized over specificity** 
for threshold selection.

### Model Behavior Over Evolutionary Time

As sequences diverge over increasingly long evolutionary timescales, the 
ROC curve is expected to degrade toward the diagonal (TPR = FPR at every 
threshold) — indicating the model becomes no better than random guessing. 
This reflects the same underlying issue as the misclassification analysis: 
extreme divergence erodes the distinguishing signal between coding and 
non-coding mutation patterns.

## Note
Datasets are simulated for coursework purposes.