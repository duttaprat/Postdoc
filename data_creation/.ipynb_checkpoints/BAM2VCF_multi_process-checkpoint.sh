#!/bin/bash

function process_bam {
    bam_file="$1"
    reference_genome="/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"

    if [ -f "$bam_file" ]; then
        echo "Processing $bam_file"
        
        # Extract the filename without the extension for the output VCF name
        base_name=$(basename "$bam_file" .bam)
        out_dir=$(dirname "$bam_file")

        # Mpileup and then call variants
        bcftools mpileup -Ou -f $reference_genome $bam_file | \
            bcftools call --ploidy 2 -mv -Oz -o $out_dir/$base_name.vcf.gz

        # Index the VCF
        bcftools index $out_dir/$base_name.vcf.gz
        
        # Delete the BAM file
        rm $bam_file
    fi
}

export -f process_bam

parent_directory="/data/projects/BAM_files/data/BAM/first"
find $parent_directory -type f -name "*.bam" | parallel -j 5 process_bam {}
