#!/usr/bin/env python
from typing import List, Tuple
import random

def upgma(
    d_matrix: List[List[int]], 
    species: List[str]
    ) -> None:
    """Performs UPGMA clustering algorithm on a distance matrix.

    Continuously finds and merges the pair with the smallest distance until only
    one cluster remains. Calls helper functions, leaving completion to students.

    Args:
        d_matrix (List[List[int]]): Symmetric distance matrix representing pairwise distances.
        species (List[str]): List of species names.

    Returns:
        None
    """
    heights = {s: 0.0 for s in species} #create a dictionary with every species height and set all to 0 to initialize


    while len(d_matrix) > 1:
        least_row: int
        least_col: int
        least_row, least_col = find_smallest(d_matrix)
        # Finds the smallest non-zero matrix coordinate
        node_height = d_matrix[least_row][least_col] / 2 #for ultrametric tree, node should be halfway between the 2 species
        branch_row = node_height - heights[species[least_row]]
        branch_col = node_height - heights[species[least_col]]

        branch_length_info = (f"({species[least_row]}:{branch_row},{species[least_col]}:{branch_col})")


        d_matrix = update_matrix(d_matrix, least_row, least_col)
        # This is the function you are writing (see below).

        species = update_species(species, least_row, least_col)
        # Updates the species list after merging a cluster.
        species[least_row] = branch_length_info
        heights[branch_length_info] = node_height

        print(d_matrix)
        print(species)

    with open("dsouza_katelyn_upgma.txt", "w") as f:
        f.write(species[0])


def find_smallest(d_matrix: List[List[int]]) -> Tuple[int, int]:
    """Finds the smallest non-zero distance in the matrix.

    Searches the matrix for the coordinate of the shortest distance between organisms.
    To complete: students should implement the logic based on assignment hints.

    Args:
        d_matrix (List[List[int]]): Symmetric distance matrix.

    Returns:
        Tuple[int, int]: Indices (row, col) of the smallest non-zero distance entry.
    """
    row: int
    col: int

    minimum_value = float('inf') #initialize a variable for the minimum value, start with infinity so any number is lower
    current_lowest = [] #create an empty list of pairs that have the mimumn value

    for i in range(len(d_matrix)): #for every row in the matrix
        for j in range(i + 1, len(d_matrix)): #for every column within that row
            value = d_matrix[i][j] #set value variable equal to the value at that coordinate point

            if value < minimum_value: #if the value is a new minimum
                minimum_value = value #replace the current_lowest list with this pair
                current_lowest = [(i,j)] #set the coordinate to the list of lowest
            elif value == minimum_value: #if the value is tied with the current minimum
                current_lowest.append((i,j)) #add the coordinate pair to the ongoing list
    row, col = random.choice(current_lowest) #choose one of the pairs in the lowest list at random and unpack the tuple
    return row, col #return the row and column of one of the lowest values


def update_species(
    species: List[str],
    row: int, 
    col: int
    ) -> List[str]:
    """Updates the species list after merging the two closest clusters.

    Combines the species labels at specified indices, deletes the merged column, and
    returns the updated list.

    Args:
        species (List[str]): List of current species/cluster names.
        row (int): Row index of cluster to merge.
        col (int): Column index of cluster to merge.

    Returns:
        List[str]: Updated species list.
    """
    species[row] = f"({species[row]},{species[col]})"
    del species[col]
    return species


def update_matrix(
    d_matrix: List[List[int]],
    row: int, 
    col: int
) -> List[List[int]]:
    """Updates the distance matrix for a new merged cluster.

    The input is the matrix and the row/col of the minimum distance pair to merge.

    Args:
        d_matrix (List[List[int]]): Previous symmetric distance matrix.
        row (int): Row index of cluster to merge.
        col (int): Column index of cluster to merge.

    Returns:
        List[List[int]]: New distance matrix with appropriate cluster distances.

    Notes:
        For more information: https://en.wikipedia.org/wiki/UPGMA
    """
    new_mat: List[List[int]] = []

    for k in range(len(d_matrix)):
        if k != row and k != col:
            average = ((d_matrix[row][k]) + (d_matrix[col][k])) / 2
            d_matrix[row][k] = average
            d_matrix[k][row] = average
    for i in range(len(d_matrix)):
        del d_matrix[i][col]
        
    del d_matrix[col]
    new_mat = d_matrix
    return new_mat


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Perform UPGMA clustering with a distance matrix and species list."
    )
    parser.add_argument(
        "--matrix",
        type=str,
        default=None,
        required=False,
        help="Path to the distance matrix file (default: built-in matrix). Format: whitespace/comma separated integers, one row per line.",
    )
    parser.add_argument(
        "--species",
        type=str,
        default=None,
        required=False,
        help="Path to the species list file (default: built-in list). Format: one species name per line.",
    )
    args = parser.parse_args()

    # Default built-in values if no files are provided
    if args.matrix:
        dist = load_distance_matrix(args.matrix)
    else:
        distance_matrix: List[List[int]] = [
            [0, 12, 12, 13, 15, 15],
            [12, 0, 2, 6, 8, 8],
            [12, 2, 0, 6, 9, 9],
            [13, 6, 6, 0, 8, 8],
            [15, 8, 9, 8, 0, 4],
            [15, 8, 9, 8, 4, 0],
        ]

    if args.species:
        sp = load_species_list(args.species)
    else:
        species_list: List[str] = [
            "M_Spacii",
            "T_Pain",
            "G_Unit",
            "Q_Doba",
            "R_Mani",
            "A_Finch",
        ]

    upgma(distance_matrix, species_list)