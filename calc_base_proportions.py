import subprocess
import glob

positions = [2466, 2473, 2485, 2521, 2527, 2535, 2562, 2597, 2615]

# Output file
output_file = "/gpfs/data/keatinglab/Twist_viral_panel_dataset/scripts/combined_output.txt"

# Open the output file in write mode and add the header
with open(output_file, "w") as f:
    f.write("BAM_File\tPosition\tA_Count\tT_Count\tC_Count\tG_Count\n")

# Iterate over all BAM files in the directory
for bam_file in glob.glob("/gpfs/data/keatinglab/Twist_viral_panel_dataset/pathseq_reports/RNA*Explant*.bam"):
    # Initialize counts for each position
    counts = {pos: {'A': 0, 'T': 0, 'C': 0, 'G': 0} for pos in positions}
    
    # Use samtools view to extract reads from the BAM file
    process = subprocess.Popen(['samtools', 'view', bam_file], stdout=subprocess.PIPE)
    
    for line in process.stdout:
        line = line.decode('utf-8')
        fields = line.split('\t')
        flag = int(fields[1])
        pos = int(fields[3])
        cigar = fields[5]
        seq = fields[9]
        
        # Skip unmapped reads
        if flag & 4:
            continue
        
        # Process CIGAR string to determine the alignment
        ref_pos = pos
        read_pos = 0
        while cigar:
            length = int(''.join(filter(str.isdigit, cigar)))
            type_ = ''.join(filter(str.isalpha, cigar))
            cigar = cigar[len(str(length)) + len(type_):]
            
            if type_ in 'M=X':
                for i in range(length):
                    if ref_pos in positions:
                        base = seq[read_pos]
                        if base in counts[ref_pos]:
                            counts[ref_pos][base] += 1
                    ref_pos += 1
                    read_pos += 1
            elif type_ in 'IS':
                read_pos += length
            elif type_ in 'DN':
                ref_pos += length
    
    # Write the results to the output file
    with open(output_file, "a") as f:
        for pos in positions:
            f.write(f"{bam_file}\t{pos}\t{counts[pos]['A']}\t{counts[pos]['T']}\t{counts[pos]['C']}\t{counts[pos]['G']}\n")

print("Processing complete.")