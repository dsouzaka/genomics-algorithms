#!/usr/bin/env python
import math
from typing import Dict, List


MODEL_CODONS: List[str] = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
    "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
    "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "CCT", "CCC",
    "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC",
    "GCA", "GCG", "TAT", "TAC", "CAT", "CAC", "CAA", "CAG",
    "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
    "TGT", "TGC", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG",
    "GGT", "GGC", "GGA", "GGG", "TGG", "TAA", "TAG", "TGA",
]
# This list is the codon labels to the coding and non-coding models.


def score_models(
    coding_probs: str = "./codingModel.tab",
    noncoding_probs: str = "./noncodingModel.tab",
    ancestor_fasta: str = "./Ancestor.fa",
    spacii_fasta: str = "./Spacii.fa"
) -> None:
    """Determines whether sequences are coding or non-coding.

    For each sequence, sets up variables and context for calculating log odds
    scores under the coding and non-coding models. Comparison and printout
    structure are present, but it is up to the student to complete the log odds
    computation.

    Args:
        coding_probs (str): Path to the probability matrix file for the coding model.
        noncoding_probs (str): Path to the probability matrix file for the non-coding model.
        ancestor_fasta (str): Path to the FASTA file containing ancestor sequences.
        spacii_fasta (str): Path to the FASTA file containing M. Spacii sequences.

    Returns:
        None
    """
    coding_matrix: List[List[float]] = get_probs(coding_probs)
    # 64x64 matrix for triplet mutations in the coding model

    noncoding_matrix: List[List[float]] = get_probs(noncoding_probs)
    # 64x64 matrix for triplet mutations in the non-coding model

    # In these matrices, the rows are ancestor codons and the columns are spacii codons

    id2ancestor_seq: Dict[str, str] = get_seq(ancestor_fasta)
    # Reads the ancestor sequences

    id2spacii_seq: Dict[str, str] = get_seq(spacii_fasta)
    # Reads the M. Spacii sequences

    all_ids: List[str] = list(id2ancestor_seq.keys())
    # The two sequence dictionaries above share indexes
    with open("dsouza_katelyn_models.txt", "w") as f: #write the results into an empty file
        for seq_id in all_ids:
            c_score: float = 0.0
            # Variable for the score using the coding model

            n_score: float = 0.0
            # Variable for the score using the non-coding model

            #######################################################################
            # YOUR CODE GOES HERE: Compute the log odds for both models.          #
            # Which sequences came from each model? In other words, which genes   #
            # are still coding for protein?                                       #
            #                                                                     #
            # Hint 1: Don't multiply probabilities and THEN take the log.         #
            # Take the log of each probability and sum. Use math.log(x).          #
            #                                                                     #
            # Hint 2: Use .index() to get the position of a codon in the list.    #
            # Example: a = ["a", "b", "c", "d"]; a.index("c") returns 2           #
            #                                                                     #
            # Hint 3: MODEL_CODONS labels row/col of matrix.                      #
            # coding_matrix[MODEL_CODONS.index("TTT")][MODEL_CODONS.index("TTA")] #
            # is the probability of TTT(ancestor) to TTA(M. Spacii) transition    #
            #######################################################################

            for i in range(0, len(id2ancestor_seq[seq_id]), 3):
                ancestor_codon = id2ancestor_seq[seq_id][i:i+3] #create variables for codons for ancestor and spacii codons
                spacii_codon = id2spacii_seq[seq_id][i:i+3]

                c_probability = coding_matrix[MODEL_CODONS.index(ancestor_codon)][MODEL_CODONS.index(spacii_codon)] #calculate the probability of coding using the coding matrix
                c_score = c_score + math.log(c_probability) #take the log of the probability
                n_probability = noncoding_matrix[MODEL_CODONS.index(ancestor_codon)][MODEL_CODONS.index(spacii_codon)] #calculate the probability of noncoding 
                n_score = n_score + math.log(n_probability) #take the log of the probability
            
            if c_score > n_score: #if there is a higher liklihood of coding than non coding
                print(
                    f"{seq_id} is coding {c_score} {n_score}" #it is coding, write to file
                )
                f.write(f"{seq_id} is coding {c_score} {n_score}\n")
            else: #otherwise it is not coding, write to file
                print(
                    f"{seq_id} is NOT coding {c_score} {n_score}"
                )
                f.write(f"{seq_id} is NOT coding {c_score} {n_score}\n")


def get_probs(filename: str) -> List[List[float]]:
    """Reads a tab-delimited probability matrix file.

    Args:
        filename (str): Path to the matrix file.

    Returns:
        List[List[float]]: Square matrix (list of lists of floats).
    """
    p_matrix: List[List[float]] = []
    with open(filename, encoding="utf-8") as file:
        for line in file:
            tmp: List[float] = [float(i) for i in line.strip().split("\t")]
            p_matrix.append(tmp)
    return p_matrix


def get_seq(filename: str) -> Dict[str, str]:
    """Reads a multi-fasta file into a dictionary.

    Args:
        filename (str): FASTA file with IDs indicated by '>'.

    Returns:
        Dict[str, str]: Dictionary mapping sequence IDs to full sequences.
    """
    id2seq: Dict[str, str] = {}
    currkey: str = ""
    with open(filename, encoding="utf-8") as file:
        for line in file:
            if line.startswith(">"):
                currkey = line[1:].split("|")[0].strip()
                id2seq[currkey] = ""
            else:
                id2seq[currkey] += line.rstrip()
    return id2seq


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Classify gene sequence pairs as coding or non-coding using "
            "probability models and aligned FASTA files."
        )
    )
    parser.add_argument(
        "--coding-probs",
        type=str,
        default="./codingModel.tab",
        help="Path to the coding model probability matrix file (default: ./codingModel.tab)",
    )
    parser.add_argument(
        "--noncoding-probs",
        type=str,
        default="./noncodingModel.tab",
        help="Path to the non-coding model probability matrix file (default: ./noncodingModel.tab)",
    )
    parser.add_argument(
        "--ancestor-fasta",
        type=str,
        default="./Ancestor.fa",
        help="Path to the ancestor sequence FASTA file (default: ./Ancestor.fa)",
    )
    parser.add_argument(
        "--spacii-fasta",
        type=str,
        default="./Spacii.fa",
        help="Path to the Spacii sequence FASTA file (default: ./Spacii.fa)",
    )

    args = parser.parse_args()
    score_models(
        coding_probs=args.coding_probs,
        noncoding_probs=args.noncoding_probs,
        ancestor_fasta=args.ancestor_fasta,
        spacii_fasta=args.spacii_fasta,
    )