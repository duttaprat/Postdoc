#!/bin/bash

parent_directory="/data/projects/BAM_files/data/BAM/first"

for subfolder in $(find $parent_directory -type d); do
    # Iterate over each BAM file in the subfolder
    for bam_file in $subfolder/*.bam; do
        if [ -f "$bam_file" ]; then
            echo "Processing $bam_file"
            
            # You can change reference genome as per your requirement
            reference_genome="/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"

            # Extract the filename without the extension for the output VCF name
            base_name=$(basename "$bam_file" .bam)
            
            # Mpileup and then call variants
            bcftools mpileup -Ou -f $reference_genome $bam_file | \
                bcftools call --ploidy 2 -mv -Oz -o $subfolder/$base_name.vcf.gz

            # Index the VCF
            bcftools index $subfolder/$base_name.vcf.gz
            
            # Delete the BAM file
            rm $bam_file
        fi
    done
done
