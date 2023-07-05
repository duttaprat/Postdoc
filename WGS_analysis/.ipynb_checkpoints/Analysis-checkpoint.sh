#!/bin/bash

# Path to your reference genome
reference="/home/pdutta/Data/Human_Genome_Data/GRCh38_latest_genomic.fna"

# Parent directory containing your BAM files and their subdirectories
parent_dir="/home/pdutta/TCGA-GDC_data/WGS/BRCA/test_data"

# Output directory for VCF files
vcf_dir="/home/pdutta/TCGA-GDC_data/WGS/BRCA/test_vcf"


# Initialize a counter for the number of BAM files
count=0


# Recursively find all BAM files in the parent directory and its subdirectories
for bam_file in $(find $parent_dir -name "*.bam")
do
  # Increment the counter
  count=$((count + 1))

  # Print the name of the BAM file and the number of files processed so far
  echo "Processing BAM file $count: $bam_file"

  # Get the base name of the file (i.e., without the directory and extension)
  base_name=$(basename $bam_file .bam)
  # Directory of bam file
  dir_name=$(dirname "$bam_file")
    
  # Subdirectory name
  sub_dir_name=$(basename "$dir_name")
  echo "BAM file is in subdirectory $sub_dir_name"

  # Record the start time
  start_time=$(date +%s)

  # Call variants with GATK and output a VCF file
  # Create a sorted version of the BAM file
  sorted_bam="$parent_dir/$sub_dir_name/${base_name}_sorted.bam"
  echo "Output Processing BAM file $count: $sorted_bam"
  samtools sort -o $sorted_bam $bam_file

  # Index the sorted BAM file
  samtools index $sorted_bam

  # Call variants and output a VCF file
  vcf_file="$vcf_dir/${base_name}.vcf"
  bcftools mpileup -Ou -f $reference $sorted_bam | bcftools call -Ov -mv - > $vcf_file
  # Record the end time
  end_time=$(date +%s)

  # Calculate and print the time taken
  time_taken=$(expr $end_time - $start_time)
  echo "Time taken for BAM file $count: $time_taken seconds \n\n"
done







