#!/usr/bin/env python3
"""
20251121_base_content.py

Minimalist HMM-based dinucleotide sequence classifier for pathogen vs. M. Spacii.
"""
from typing import Dict, List
import math
import matplotlib.pyplot as plt


BASE_IDX: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "T": 3}


def get_seq(
    fasta_filename: str,
) -> Dict[str, str]:
    """
    Load sequences from a FASTA file into a dictionary.

    Args:
        fasta_filename (str): Path to a FASTA-format file.

    Returns:
        Dict[str, str]: Dictionary mapping sequence IDs to nucleotide strings.
    """
    seqs: Dict[str, str] = {}
    current_id: str = ""
    with open(fasta_filename, "r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                current_id = line[1:]
                seqs[current_id] = ""
            else:
                seqs[current_id] += line
    return seqs


def train_model(
    model: List[List[float]],
    fasta_filename: str,
) -> List[List[float]]:
    """
    Train a dinucleotide Markov model on the sequences in the FASTA file.
    
    This should look at all the training data and calculate how many
    dinucleotides preceed each base similar to what was outlined on the slides

    Args:
        model (List[List[float]]): Preallocated 4x4 matrix.
        fasta_filename (str): Path to the training FASTA file.

    Returns:
        model (List[List[float]]): The normalized transition probability matrix.
    """
    seqs = get_seq(fasta_filename)

    for i in range(4):
        for j in range(4):
            model[i][j] = 1.0

    for seq in seqs.values(): #for every sequence
        seq = seq.upper() #make sure the sequence is all caps
        for k in range(len(seq) - 1): #for every base in the sequence stopping 1 before the end
            b1 = seq[k] #set b1 equal to the current base
            b2 = seq[k+1] #set b2 equal to the next base 
            if b1 in BASE_IDX and b2 in BASE_IDX: #if the base in b1 is a proper base
                i = BASE_IDX[b1] #convert letter into number 0-3
                j = BASE_IDX[b2] #convert letter into number 0-3
                model[i][j] += 1.0 #every time you see the nucleotide in training sequence add 1
    for i in range(4): #for each row with raw counts
        row_sum = sum(model[i]) #add all the counts in the row and divide every cell by it so each row sums to 1 and represents the probability of the next base given the current base
        for j in range(4):
            model[i][j] /= row_sum
        
    return model 


def get_log_likelihood(
    model1: List[List[float]],
    model2: List[List[float]],
    seq: str,
) -> float:
    """
    Compute the log-likelihood ratio for the sequence under two Markov models.

    Args:
        model1 (List[List[float]]): First trained 4x4 model (numerator).
        model2 (List[List[float]]): Second trained 4x4 model (denominator).
        seq (str): Nucleotide sequence.

    Returns:
        score (float): log(P(seq | model1)) - log(P(seq | model2)).
    """
    score = 0.0 #initialize the running total to start at 0
    seq = seq.upper() #make sure all bases are uppercase

    for k in range(len(seq) - 1): #for every position in the sequence stopping 1 before the end
        b1 = seq[k]
        b2 = seq[k + 1] #take the current base and the one after
        if b1 in BASE_IDX and b2 in BASE_IDX: #if both are real bases
            i = BASE_IDX[b1] #convert the bases to numbers for matrix positions
            j = BASE_IDX[b2]
            score += math.log(model1[i][j]) - math.log(model2[i][j]) #to find out how probable the pair is in both models take the log of each and add the difference to score, if Spacii model1 asigns higher probabilities than pathogen model2 the score is positive, if model 2 is more then score is negative
    return score



def scores_output_text(
    spacii_scores: List[float],
    pathogen_scores: List[float],
    spacii_ids: List[str],
    pathogen_ids: List[str],
    output_filename: str = "results.tab",
) -> None:
    """
    Write Spacii and pathogen log-likelihood ratios to a tab-delimited file.

    Args:
        spacii_scores (List[float]): Scores for Spacii test sequences.
        pathogen_scores (List[float]): Scores for pathogen test sequences.
        output_filename (str): Output filename.

    Returns:
        None
    """
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write("SeqID\tClassification\tScore\n")
        for seq_id, score in zip(spacii_ids, spacii_scores):
            f.write(f"{seq_id}\tM. Spacii\t{score}\n")
        for seq_id, score in zip(pathogen_ids, pathogen_scores):
            f.write(f"{seq_id}\tpathogen\t{score}\n")


def main(
    spacii_fa: str,
    pathogen_fa: str,
    spacii_fa_t: str,
    pathogen_fa_t: str,
    results_tab: str = "results.tab",
    histogram_png: str = "hmm_score_histogram.png",
) -> None:
    """
    High-level workflow: load data, train models, score sequences, plot, and save results.

    Args:
        spacii_fa (str): Path to M. Spacii test FASTA.
        pathogen_fa (str): Path to pathogen test FASTA.
        spacii_fa_t (str): Path to M. Spacii training FASTA.
        pathogen_fa_t (str): Path to pathogen training FASTA.
        results_tab (str): Output tab file for results.
        histogram_png (str): Output PNG file for histogram.

    Returns:
        None
    """
    spacii_id2seq = get_seq(spacii_fa)
    pathogen_id2seq = get_seq(pathogen_fa)

    spacii_train_model = [[0.0 for _ in range(4)] for _ in range(4)]
    path_train_model = [[0.0 for _ in range(4)] for _ in range(4)]

    spacii_train_model = train_model(spacii_train_model, spacii_fa_t)
    path_train_model = train_model(path_train_model, pathogen_fa_t)

    markov_scores_spacii = [
        get_log_likelihood(spacii_train_model, path_train_model, seq)
        for seq in spacii_id2seq.values()
    ]
    markov_scores_path = [
        get_log_likelihood(spacii_train_model, path_train_model, seq)
        for seq in pathogen_id2seq.values()
    ]

    plt.hist(
        [markov_scores_path, markov_scores_spacii],
        bins=20,
        label=["pathogen", "spacii"],
        rwidth=1.0,
        density=True,
    )
    plt.legend()
    plt.xlabel("Log-likelihood ratio (Spacii - Pathogen)")
    plt.ylabel("Density")
    plt.savefig(histogram_png)
    scores_output_text(
        markov_scores_spacii, markov_scores_path,
        list(spacii_id2seq.keys()), list(pathogen_id2seq.keys()),
        results_tab
    )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="HMM dinucleotide classifier: pathogen vs. M. Spacii"
    )
    parser.add_argument(
        "--spacii-fa",
        type=str,
        default="MSpacii.fa",
        help="FASTA file with M. Spacii sequences (test)",
    )
    parser.add_argument(
        "--pathogen-fa",
        type=str,
        default="pathogen.fa",
        help="FASTA file with pathogen sequences (test)",
    )
    parser.add_argument(
        "--spacii-fa-t",
        type=str,
        default="MSpacii_training.fa",
        help="FASTA file with M. Spacii training sequences",
    )
    parser.add_argument(
        "--pathogen-fa-t",
        type=str,
        default="pathogen_training.fa",
        help="FASTA file with pathogen training sequences",
    )
    parser.add_argument(
        "--results-tab",
        type=str,
        default="results.tab",
        help="Output tab file for classification results",
    )
    parser.add_argument(
        "--histogram-png",
        type=str,
        default="hmm_score_histogram.png",
        help="Output PNG for score histogram",
    )

    args = parser.parse_args()
    main(
        spacii_fa=args.spacii_fa,
        pathogen_fa=args.pathogen_fa,
        spacii_fa_t=args.spacii_fa_t,
        pathogen_fa_t=args.pathogen_fa_t,
        results_tab=args.results_tab,
        histogram_png=args.histogram_png,
    )