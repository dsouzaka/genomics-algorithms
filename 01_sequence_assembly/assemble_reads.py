def read_sequences(filename):
    """read sequences from a text file (one sequence per line)"""
    sequences = [] #create an empty list to hold the sequences read from file
    with open(filename, 'r') as f: #opening the file in read mode as a file object f and automatically closes file when done
        for line in f: #loop through each line in f, for each line
            seq = line.strip() #remove the whitespace from the ends
            if seq: #if the line exists 
                sequences.append(seq) #add the sequence to the end of the list of sequences
    return sequences #return the list of sequences

def find_overlaps (a, b, k):
    """find the overlap between sequences a and b with k as the mimimum overlap required"""
    max_overlap = min(len(a), len(b)) #can't overlap more than the smaller of the 2 
    for overlap_len in range(max_overlap, k-1, -1): #for each possible overlap length starting with the max and ending with the mimimum 10 by default
        if a[-overlap_len:] == b[:overlap_len]: #for each length x, if the last x characters of a match the first x characters of b
            return overlap_len #return that overlap length, keep looping through until it does to get the largest possible overlap
    return 0 #if it never overlaps, return 0; also returns 0 for anything less than k

def write_contigs(contigs, output_file):
    """write assembled contigs to FASTA file"""
    with open(output_file, 'w') as f: #open file in write mode
        for i, contig in enumerate(contigs): #for every contig and its index
            f.write(f">Contig{i+1}\n")  #create a FASTA header that is the index+1 so the first contig is 1 and not 0
            f.write(f"{contig}\n") #write the sequence below 

def assemble_reads(inputfile, outputfile, k=10):
    """assemble reads into contigs"""
    reads = read_sequences(inputfile) #read all sequences from file
    print(f"Loaded {len(reads)} reads") #tells you how many reads to ensure it was read properly
    used = [False] * len(reads) #create a list of boolean values that correspond to the list of reads to track which have been incorporated, start as FALSE by default until incorporated
    contigs = [] #create an empty list to store created contigs
    
    while not all(used): #as long as all the values in the "used" list are not true
        contig = None #create a contig variable and set to None since we don't have one yet
        for i in range(len(reads)): #for every read 
            if not used[i]: #if it has not been used yet
                contig = reads[i] #start a contig with this read
                used[i] = True #mark it as used
                break #stop since we have a read
        
        extended = True #keeps track of if extension can keep occuring
        while extended: #while extended is true
            extended = False #set to false until we find overlap
            best_overlap = 0 #set the best overlap to 0 until wefind one
            best_idx = None #set the best overlap id to None until we find one
            extend_direction = None #tracks if we are extending to the right or left, None until we start extending

            for i in range(len(reads)): #for every unused read
                if not used[i]: #for every unused read
                    # first try to extend right
                    overlap_right = find_overlaps(contig, reads[i], k) #checks if the read overlaps by at least k to the right of the current contig
                    if overlap_right > best_overlap: #if it overlaps more than the best overlap
                        best_overlap = overlap_right #set the overlap as the new best overlap
                        best_idx = i #set the best id to the id
                        extend_direction = 'right' #set extend direction to right
                    
                    # then try to extend left
                    overlap_left = find_overlaps(reads[i], contig, k) #check is the read overlaps by at least k to the left of the current contig
                    if overlap_left > best_overlap: #if it overlaps more than the best overlap
                        best_overlap = overlap_left #set the overlap as the new best overlap
                        best_idx = i #set the best id to the id
                        extend_direction = 'left' #set extend direction to left
            
            if best_idx is not None: #if we found an overlap
                if extend_direction == 'right': #if this overlap is to the right
                    contig += reads[best_idx][best_overlap:] #add the non-overlapping part to the right (ensures we dont repeat sections)
                else: #otherwise extend to the left
                    contig = reads[best_idx][:-best_overlap] + contig #add the non-overlapping part to the left
                
                used[best_idx] = True #set this read as used
                extended = True #set extended to true, keep looping until we cant extend anymore
        
        contigs.append(contig) #after extension is done, add contig to the list of contigs
        print(f"Assembled contig {len(contigs)}: length {len(contig)} bp") #prints the number and length of each contig as it is created
    
    write_contigs(contigs, outputfile) #write the contigs to FASTA file
    return contigs, reads #function returns contigs and reads

def calculate_read_depth(contigs, reads):
    """calculate read depth for each position in each contig"""
    contig_depths = [] #create an empty list for contig read depths
    
    for contig_idx, contig in enumerate(contigs): #for each contig
        depth = [0] * len(contig) #create an array of zeros for each position to record depth
        
        for read in reads: #for every read
            for pos in range(len(contig) - len(read) + 1): #try every position and check if the read matches
                if contig[pos:pos+len(read)] == read: #if it matches read
                    for i in range(len(read)): #add 1 to all positions
                        depth[pos + i] += 1
                    break #move to next read
        
        contig_depths.append(depth) #save the contigs depth array
    return contig_depths

def calculate_density(depth_list):
    """calculate density from depth counts"""
    total_depth = sum(depth_list) #add all the depth values for each contig
    if total_depth == 0: #if the total depth is 0, return 0
        return [0] * len(depth_list)
    return [d / total_depth for d in depth_list] #otherwise divide each depth by total depth to get density

def export_to_csv(contigs, contig_depths, contig_number=0, output_file="read_density_data.csv"):
    """export read depth and density data to CSV for excel"""
    import csv
    
    with open(output_file, 'w', newline='') as f: #open output file in write mode
        writer = csv.writer(f) 
        
        writer.writerow(['Position', f'Contig{contig_number+1}', f'density {contig_number+1}']) #create the header that shows position, contig number, and density for contig number
        
        density = calculate_density(contig_depths[contig_number]) #calculate the density for the contig
        
        for pos in range(len(contig_depths[contig_number])): #for each position's data, write the data into csv
            depth = contig_depths[contig_number][pos]
            dens = density[pos]
            writer.writerow([pos+1, depth, f'{dens:.6f}']) #make sure position starts at 1
        
        writer.writerow(['', '', '']) #add a blank row before the sum
        writer.writerow(['sum', sum(contig_depths[contig_number]), ''])
    
    print(f"\nData for Contig {contig_number+1} exported to {output_file}") #print the name of the output file for csv



contigs, reads = assemble_reads("seqReadFile2023.txt", "dsouza_katelyn_ps1.fasta", k=10) #use assemble reads function to run the assembly 

contig_depths = calculate_read_depth(contigs, reads) #use the calculate read depth function to calculate read depth

export_to_csv(contigs, contig_depths, contig_number=0, output_file="read_density_data.csv") #export to csv for excel