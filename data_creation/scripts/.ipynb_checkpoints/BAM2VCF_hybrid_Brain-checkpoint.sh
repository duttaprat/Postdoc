chmod +x /home/campus.stonybrook.edu/pdutta/gdc-client
export PATH=$PATH:/home/campus.stonybrook.edu/pdutta




export CANCER_TYPE="Brain" 
export REF_PATH="/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"
export MANIFEST_PATH="/data/projects/GDC_Cancer_Wise/Brain/Manifest/Brain_gdc_manifest.2023-08-30_2.txt"
export TOKEN_PATH='/data/projects/BAM_files/token/gdc-user-token.2023-08-18T17_34_29.807Z.txt'
export DOWNLOAD_PATH="/data/projects/GDC_Cancer_Wise/$CANCER_TYPE/Data/BAM_multi/"
export VCF_PATH="/data/projects/GDC_Cancer_Wise/$CANCER_TYPE/Data/VCF_multi/"

chmod 600 $TOKEN_PATH

cd ..
python GDC_client_BAM2VCF_hybrid.py $TOKEN_PATH $REF_PATH $MANIFEST_PATH $DOWNLOAD_PATH $VCF_PATH