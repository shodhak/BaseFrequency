#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --job-name=nucleotide_proportions
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --nodes=1

module load samtools
module load bedtools

# Define the reference genome and positions of interest
REFERENCE_GENOME="/gpfs/data/keatinglab/Twist_viral_panel_dataset/appv_genome/appv_ExplantKidney/appv_ExplantKidney.fasta"

cd /gpfs/data/keatinglab/Twist_viral_panel_dataset/pathseq_reports

# Create or overwrite the output TSV file and add the header
output_file="/gpfs/data/keatinglab/Twist_viral_panel_dataset/scripts/coverage_depth_report.tsv"
echo -e "BAM_File\tAverage_Depth\tPercent_Genome_Coverage" > $output_file

for bam_file in *_sorted_Bowtie2_ExplantKidney.bam
do
    file_name=${bam_file%%.*}
    echo "Processing $file_name"
    
    # Calculate average depth
    average_depth=$(samtools depth $bam_file | awk '{sum+=$3} END { print sum/NR}')
    echo "Average Depth = $average_depth"
    
    # Calculate percent genome coverage
    genome_coverage=$(bedtools genomecov -ibam $bam_file | awk '{if($1=="genome" && $2=="0") uncovered+=$3; else covered+=$3} END {print (covered/(covered+uncovered))*100}')
    echo "Percent Genome Coverage = $genome_coverage%"
    
    # Append the results to the output file
    echo -e "${bam_file}\t${average_depth}\t${genome_coverage}%" >> $output_file
done