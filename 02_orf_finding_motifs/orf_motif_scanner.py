from typing import Optional, List, Tuple
import random
import csv

BASE2IDX: dict[str, int] = {base: idx for idx, base in enumerate("ATCG")}

MOTIF: list[list[float]] = [
    [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0,   2, -99, -99, 0.5], # A
    [  0,   0,   0,   0, 0, 0, 0, 0, 0, -99,   2, -99,   0], # T
    [  0,   0,   0,   0, 0, 0, 0, 0, 0, -99, -99, -99,   0], # C
    [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0.5, -99,   2,   0]  # G
]

def scan_seq(
    sequence: str,
    start_codons: Optional[List[str]] = None,
    stop_codons: Optional[List[str]] = None,
    min_length: int = 60,
) -> List[Tuple[str, int, int]]:
    """Find potential open reading frames (ORFs) in the given sequence and return positional information.

    Args:
        sequence (str): Sequence of contiguous nucleotides.
        start_codons (Optional[List[str]]): List of canonical or non-canonical start codons. Defaults to None.
        stop_codons (Optional[List[str]]): List of canonical or non-canonical stop codons. Defaults to None.
        min_length (int): Minimum length of the ORF. Defaults to 60.

    Returns:
        List[Tuple[str, int, int]]: Tuples containing the ORF sequence, start position relative to the contig, and ORF length.
    """
    #set defauly start and stop codons to use if not provided
    if start_codons is None:
        start_codons = ["ATG", "GTG"]
    if stop_codons is None:
        stop_codons = ["TAA", "TAG", "TGA"]
    
    #create an empty list to store found ORFs
    orfs = []

    for frame in range(3): #scan all 3 reading frames
        i = frame #set i to the start position for the frame being checked
        while i <= len(sequence) - 3: #need at least 3 bases for a codon
            codon = sequence[i:i+3] #set codon to the 3 bases

            if codon in start_codons: #if codon is a start codon, set index as start position
                start_pos = i

                j = i+3 #start searching for stop codon after the start
                while j <= len(sequence) - 3: #for every codon after, set as potential stop codon
                    stop_codon = sequence[j:j+3]

                    if stop_codon in stop_codons: #if it actually is a stop codon, this is an ORF
                        orf_seq = sequence[start_pos : j+3] #extract the sequence and save to to orf_seq variable including the stop codon
                        orf_length = len(orf_seq) #save the length of the sequence to orf_length

                        if orf_length >= min_length and orf_length % 3 == 0: #if it meets the minimum length and is divisible by 3 
                            orfs.append((orf_seq, start_pos, orf_length)) #add tuple to orfs list

                        i = j + 3 #after finding ORF, go past it and continue scanning
                        break
                    j += 3
                else: #move forward by 3 if inner loop does not find stop codon
                    i += 3
                    break
            else: #if no start codon found, move to the next codon
                i += 3
    return orfs


def score_motif(
    motif_seq: str,
    motif_scoring_matrix: Optional[List[List[float]]] = None,
) -> float:
    """Generate a motif score for the provided sequence using the given scoring matrix.

    Args:
        motif_seq (str): 13 base pair motif sequence.
        motif_scoring_matrix (Optional[List[List[float]]]): 2D position-by-base score matrix.
            If None, uses global MOTIF.

    Returns:
        float: Motif score calculated for the sequence.
    """
    #default to MOTIF if scoring matrix not provided
    if motif_scoring_matrix is None:
        motif_scoring_matrix = MOTIF
    
    #initialize the score to 0
    score = 0.0

    #score each position 0-12
    for position, base in enumerate(motif_seq): #convert each base letter to a row number for matrix 
        row = BASE2IDX[base] #get row index, A=0, T=1, C=2, G=3
        score += motif_scoring_matrix[row][position] #look up the score in the matrix and add it to the total
    return score


def identify_ORFs(
    fasta_input: str, #name of file input
    score_cutoff: float, #minimum score ORF needs (default 7.25)
    min_length: int, #minimum ORF length (default 60)
    output_filename: str, #name of file for results
    start_codons: Optional[List[str]] = None, #custom start codons (defaults to ATG/GTG)
    stop_codons: Optional[List[str]] = None, #custom stop codons (defaults to TAG, TAA, TGA)
) -> None:
    """Collect candidate ORFs from contigs and write them to a specified file.

    Args:
        fasta_input (str): FASTA file containing contigs to analyze.
        score_cutoff (float): Threshold score to determine if a candidate ORF is selected.
        min_length (int): Minimum allowed length for ORF candidates.
        output_filename (str): Filename to write selected ORFs in FASTA format.
        start_codons (Optional[List[str]]): Custom start codons. Defaults to None.
        stop_codons (Optional[List[str]]): Custom stop codons. Defaults to None.

    Returns:
        None
    """
    if start_codons is None:
        start_codons = ["ATG", "GTG"]
    if stop_codons is None:
        stop_codons = ["TAA", "TAG", "TGA"]

    contigs = {} #create a dictionary with contig name and sequence

    with open(fasta_input, 'r') as f: #open fasta file in read mode
        current_header = None #initialize variable to store name of contig
        current_seq = [] #create empty list to store sequence for current contig

        for line in f: #for every line in the file
            line = line.strip() #remove the whitespace
            if line.startswith('>'): #if line starts with > it is a header
                if current_header: #save previous contig, first time gets skipped
                    contigs[current_header] = ''.join(current_seq) #join all the sequence lines into one string for contig
                current_header = line[1:].strip() #remove >
                current_seq = [] #empty the current seq variable for new contig
            else: #if its not a header line, add it to the list for the sequence
                current_seq.append(line)

        if current_header: #last contig in file gets saved to dictionary since there is no new header after it to initialize the save
            contigs[current_header] = ''.join(current_seq)
    
    with open(output_filename, 'w') as out: #write results to file
        for contig_name, sequence in contigs.items(): 
            # Use scan_seq to find all ORFs in contig
            orfs = scan_seq(sequence, start_codons, stop_codons, min_length)
            
            for orf_seq, start_pos, orf_length in orfs:
                # extract 13bp context, 9 bases before start, start codon, and 1 base after
                if start_pos >= 9 and start_pos + 4 <= len(sequence):
                    context_seq = sequence[start_pos - 9 : start_pos + 4]
                    
                    if len(context_seq) == 13: #score the context so that start codon is positions 9-10-11
                        motif_score = score_motif(context_seq)
                        
                        if motif_score >= score_cutoff:
                            #write the ORF to output including its name, length, and start position
                            header = f"> contig {contig_name} | Length {orf_length} | at position {start_pos}\n"
                            out.write(header)
                            out.write(orf_seq + "\n")


def generate_reads(
    read_length: int = 20,
    n_reads: int = 20,
    alphabet: str = "ACGT",
) -> List[str]:
    """Generate a list of synthetic reads from the given alphabet.

    Args:
        read_length (int): Length of each read to be generated. Defaults to 20.
        n_reads (int): Number of reads to generate. Defaults to 20.
        alphabet (str): String indicating allowed nucleotide bases. Defaults to "ACGT".

    Returns:
        List[str]: List of random DNA/RNA sequences.
    """
    return ["".join(random.choices(alphabet, k=read_length)) for _ in range(n_reads)]

def check_multiple_orfs(fasta_input, min_length=60): #find sequences with multiple ORFs
    contigs = {}
    with open(fasta_input, 'r') as f: #read FASTA file and add sequences to contigs dictionary
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    contigs[current_header] = ''.join(current_seq)
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            contigs[current_header] = ''.join(current_seq)
    
    print("Sequences with Multiple ORFs:\n") #find sequences with multiple ORFs
    for contig_name, sequence in contigs.items(): #loop through the contig dictioary 
        orfs = scan_seq(sequence, min_length=min_length) #find ORFs in sequence using scan_seq() function
        if len(orfs) > 1: #if the sequence has more than 1 ORF
            print(f"{contig_name}: {len(orfs)} ORFs") #print the contig name and the number of ORFs
            for i, (orf_seq, start_pos, orf_length) in enumerate(orfs, 1): #loop through each ORF and gets each tuple
                frame = start_pos % 3 #calculate which reading frame the ORF is in
                
                # extract and score context
                if start_pos >= 9 and start_pos + 4 <= len(sequence): #make sure there is enough space to extract the 13 bp 
                    context_seq = sequence[start_pos - 9 : start_pos + 4] #extract the 13bp around start codon so that start codon is positions 9 10 and 11
                    if len(context_seq) == 13: #if the length is 13, score it
                        score = score_motif(context_seq)
                    else: #if not set score to 0
                        score = 0.0
                else: #if the start position is not in the right place to extract context, set to 0
                    score = 0.0
                    
                print(f"  ORF{i}: pos {start_pos}, length {orf_length}, frame {frame}, score {score:.1f}") #print details for each ORF
            print()

def analyze_borderline_orfs(fasta_input, score_cutoff=7.25, min_length=60):
    """Find ORFs that scored close to 7.25 (6.5-7.24) and almost passed filters"""
    contigs = {}
    with open(fasta_input, 'r') as f: #read FASTA
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    contigs[current_header] = ''.join(current_seq)
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            contigs[current_header] = ''.join(current_seq)
    
    print("ORFs Close to Cutoff (scored 6.5-7.24):\n")
    
    borderline = [] #create an empty list for the borderline scoring ORFs
    
    for contig_name, sequence in contigs.items(): #for each sequence find all ORFs using scan_sequence() function
        orfs = scan_seq(sequence, min_length=min_length)
        
        for orf_seq, start_pos, orf_length in orfs:
            #extract and score context
            if start_pos >= 9 and start_pos + 4 <= len(sequence):
                context_seq = sequence[start_pos - 9 : start_pos + 4]
                
                if len(context_seq) == 13:
                    score = score_motif(context_seq)
                    
                    if 6.5 <= score < score_cutoff: #if the score is in the range to be borderline
                        borderline.append((contig_name, start_pos, orf_length, score, orf_seq)) #save the ORF as a tuple in the borderline list, containing its name, position, length, score, and the actual sequence
                        print(f"{contig_name}:") #print details for borderline ORFs for analysis 
                        print(f"  Position: {start_pos}")
                        print(f"  Length: {orf_length} bp")
                        print(f"  Score: {score:.2f} (missed by {score_cutoff - score:.2f})")
                        print(f"  Context 13bp: {context_seq}")
                        print(f"  Starts: {orf_seq[:3]}, Ends: {orf_seq[-3:]}")
                        print()
    
    if not borderline: #if there are no borderline ORFs, print that
        print("No ORFs close to the cutoff found.")
    else: #if there are, print the number of borderline ORFs found total
        print(f"Found {len(borderline)} borderline ORF(s)")
    
    return borderline #return the list with tuples

def find_orfs_in_random(sequence):
    """ Find all ORFs in a random sequence
    ATG start, TAA, TAG, TGA as stop
    """
    start_codon = "ATG" #only look for ATG as start codon
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = [] #create empty list for ORF lengths

    for frame in range(3): #for all the reading frames
        i = frame
        while i <= len(sequence) - 3: #loop through 3 bases at a time and extract each codon
            codon = sequence[i:i+3]

            if codon == start_codon: #if an ATG is found, mark its position
                start_pos = i
                j = i + 3 

                while j <= len(sequence) - 3: #keep scanning sequence to search for stop codon
                    stop_codon = sequence[j:j+3]

                    if stop_codon in stop_codons: #if a stop codon is found
                        orf_length = j + 3 - start_pos #calculate ORF length and save it 
                        orfs.append(orf_length)
                        i = j + 3 
                        break
                    j += 3 #if not a stop codon, go to next in frame
                else:
                    i += 3 #if no stop is found, break
                    break
            else: #if not a start codon, move to next codon
                i += 3
    return orfs #return list of ORF lengths

def part3_analysis(): #generate 1000 random sequences, find ORFs, analyze distribution, export the data
    n_sequences = 1000 #generate 1000 sequences
    read_length = 150 #each 150 bp long

    all_orf_lengths = [] #create variable to hold ORF lengths for all sequences
    sequences_with_long_orfs = 0 #count sequences that have ORFs longer than 60 bp

    for n in range(n_sequences): #for each sequence
        seq = generate_reads(read_length = read_length, n_reads = 1)[0] #generate 1 random 150 bp sequence
        orfs = find_orfs_in_random(seq) #find the ORFs in this sequence
        all_orf_lengths.extend(orfs) #add ORF lengths from sequence to list

        if any(length > 60 for length in orfs): #if the sequence has ORF > 60 bp, add to counter
            sequences_with_long_orfs += 1
    print(f"Generating Random Sequences and Analysis") #print stats found
    print(f"Total sequences generated: {n_sequences}")
    print(f"Total ORFs found: {len(all_orf_lengths)}")

    if all_orf_lengths:
        print(f"Mean ORF length: {sum(all_orf_lengths)/len(all_orf_lengths):.1f} bp")
    
    fraction = sequences_with_long_orfs / n_sequences
    print(f"\nSequences with ORF >60bp: {sequences_with_long_orfs}/{n_sequences}")
    print(f"Fraction: {fraction:.4f} ({fraction*100:.2f}%)")

    with open('orf_data.csv', 'w', newline='') as f: #export ORF lengths to CSV file
        writer = csv.writer(f)
        writer.writerow(['ORF_Length_bp'])
        for length in all_orf_lengths:
            writer.writerow([length])
    
    print(f"\n Data exported to 'orf_data.csv'")


if __name__ == "__main__": #run code when file is run directly
    identify_ORFs(
        fasta_input="spaceSeq.fa", 
        score_cutoff=7.25,
        min_length=60,
        output_filename="dsouza_katelyn_ps2.fa"
    ) #run main analysis to read input file, find ORFs that meet criteria, and write results to output file
    print("Output written to dsouza_katelyn_ps2.fa")

    check_multiple_orfs("spaceSeq.fa", min_length=60) #run analysis to find sequences with multiple ORFs
    borderline = analyze_borderline_orfs("spaceSeq.fa", score_cutoff=7.25, min_length=60) #run analysis to find borderline ORFs
    part3_analysis() #run the simulation to generate random sequences and find ORFs
 