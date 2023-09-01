export CANCER_TYPE="Brain" 
export REF_PATH="/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"
export BAM_PATH="/data/projects/GDC_Cancer_Wise/"$CANCER_TYPE"/Data/BAM_NEED_TO_BE_PROCESSED/"
#export BAM_PATH="/data/projects/GDC_Cancer_Wise/Brain/Data/BAM_multi/b8ee1244-f7ed-4816-953d-f1c55efba0d8"
export VCF_PATH="/data/projects/GDC_Cancer_Wise/"$CANCER_TYPE"/Data/VCF_multi/"


cd ..
python BAM2VCF_chrom.py $REF_PATH $BAM_PATH $VCF_PATH