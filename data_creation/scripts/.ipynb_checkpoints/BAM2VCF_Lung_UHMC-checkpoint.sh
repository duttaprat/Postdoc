export CANCER_TYPE="Lung" 
export REF_PATH="/home/pdutta/Rekha_LabWork/Collaborator_work/RNA-seq/RNAseq_CommonFiles/STAR/Genome/GRCh38.primary_assembly.genome.fa"
export BAM_PATH="/home/pdutta/Data/Cancer_wiseGDC/Data/"$CANCER_TYPE"/NEW_BAM/12/"
#export BAM_PATH="/data/projects/GDC_Cancer_Wise/Brain/Data/BAM_multi/b8ee1244-f7ed-4816-953d-f1c55efba0d8"
export VCF_PATH="/home/pdutta/Data/Cancer_wiseGDC/Data/"$CANCER_TYPE"/VCF/"


cd ..
python BAM2VCF_chrom.py $REF_PATH $BAM_PATH $VCF_PATH