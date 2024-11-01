#!/bin/bash

module load samtools

# Function to print error and exit
error_exit() {
    echo "ERROR: $1" >&2
    exit 1
}

# Check if required tools are installed
for tool in samtools awk; do
    if ! command -v "$tool" &> /dev/null; then
        error_exit "$tool is not installed. Please install it first."
    fi
done

# Check number of arguments
if [ $# -ne 1 ]; then
    error_exit "Usage: $0 <input.bam>"
fi

# Input file
BAM_FILE="$1"

# Validate input file exists and is readable
[ -f "$BAM_FILE" ] || error_exit "BAM file does not exist: $BAM_FILE"
[ -r "$BAM_FILE" ] || error_exit "BAM file is not readable: $BAM_FILE"

# Validate BAM file is indexed
if [ ! -f "${BAM_FILE}.bai" ]; then
    echo "BAM file is not indexed. Attempting to index..."
    samtools index "$BAM_FILE" || error_exit "Failed to index BAM file"
fi

# Array of positions to analyze
positions=("2466" "2473" "2485" "2521" "2527" "2535" "2562" "2597" "2615")

# Create output and debug files
OUTPUT_FILE="${BAM_FILE/_pathseq_output_reads_sorted_Bowtie2_ExplantKidney.bam/_base_freq.txt}"
DEBUG_FILE="base_frequencies_debug.txt"
> "$OUTPUT_FILE"
> "$DEBUG_FILE"

# Header for output file
echo "Position | A | C | G | T | Total_Depth" >> "$OUTPUT_FILE"
echo "Detailed Debugging Information" >> "$DEBUG_FILE"

# Flag to track if any data was found
data_found=0

# Check if BAM contains reads
total_reads=$(samtools view -c "$BAM_FILE")
echo "Total reads in BAM file: $total_reads" >> "$DEBUG_FILE"
if [ "$total_reads" -eq 0 ]; then
    error_exit "No reads found in the BAM file."
fi

# Analyze each position
for pos in "${positions[@]}"; do
    echo "Analyzing position $pos" >> "$DEBUG_FILE"
    
    # Use samtools view to extract reads at the specific position
    reads=$(samtools view "$BAM_FILE" | awk -v POS="$pos" '
        BEGIN {
            a_count = 0
            c_count = 0
            g_count = 0
            t_count = 0
        }
        {
            # Split CIGAR string to understand read alignment
            cigar = $6
            split($10, seq, "")  # Split sequence into individual bases
            
            # Track read position and reference position
            read_pos = 1
            ref_pos = int($4)
            
            # Parse CIGAR to find exact base position
            split(cigar, cigar_ops, /[0-9]+/)
            split(cigar, cigar_nums, /[MIDNSHP=X]/)
            
            for (i = 1; i < length(cigar_nums); i++) {
                num = cigar_nums[i]
                op = cigar_ops[i+1]
                
                # Check if position is within this read
                if (op == "M") {  # Matches/mismatches
                    if (ref_pos + num >= POS && ref_pos <= POS) {
                        # Calculate offset
                        offset = POS - ref_pos
                        if (offset >= 0 && offset < length($10)) {
                            base = tolower(seq[offset+1])
                            if (base == "a") a_count++
                            else if (base == "c") c_count++
                            else if (base == "g") g_count++
                            else if (base == "t") t_count++
                        }
                    }
                    ref_pos += num
                    read_pos += num
                }
                # Handle other CIGAR operations (insertions, deletions, etc.)
                else if (op == "I") {
                    read_pos += num
                }
                else if (op == "D") {
                    ref_pos += num
                }
            }
        }
        END {
            printf "%s|%d|%d|%d|%d\n", POS, a_count, c_count, g_count, t_count
        }
    ')
    
    # Get total depth for the position
    depth=$(samtools depth -r "$pos" "$BAM_FILE" | awk '{print $3}')
    
    # Combine reads and depth
    full_info=$(echo "$reads" | awk -v depth="$depth" -F'|' '{printf "%s|%s|%d\n", $1, $0, depth}')
    
    # Debug output
    echo "Reads for position $pos: $reads" >> "$DEBUG_FILE"
    echo "Depth for position $pos: $depth" >> "$DEBUG_FILE"
    
    # Check if we found any data
    if [ -n "$full_info" ]; then
        echo "$full_info" >> "$OUTPUT_FILE"
        data_found=1
    else
        echo "No data found for position $pos" >&2
        echo "No data found for position $pos" >> "$DEBUG_FILE"
    fi
done

# Final check for data
if [ $data_found -eq 0 ]; then
    error_exit "No data found for any of the specified positions. Check your BAM file and positions."
fi

echo "Base frequency analysis complete. Results in $OUTPUT_FILE"
echo "Detailed debug information in $DEBUG_FILE"