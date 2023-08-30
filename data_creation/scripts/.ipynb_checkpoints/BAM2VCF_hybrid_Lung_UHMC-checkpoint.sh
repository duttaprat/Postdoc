chmod +x /home/pdutta/gdc-client
export PATH=$PATH:/home/pdutta




export CANCER_TYPE="Lung" 
export REF_PATH="/home/pdutta/Rekha_LabWork/Collaborator_work/RNA-seq/RNAseq_CommonFiles/STAR/Genome/GRCh38.primary_assembly.genome.fa"
export MANIFEST_PATH="/home/pdutta/Data/Cancer_wiseGDC/Manifest_files/BronchusLung_gdc_manifest.2023-08-30_2.txt"
export TOKEN_PATH='/home/pdutta/Data/GDC/RAMANA/Token/gdc-user-token.2023-08-18T17_34_29.807Z.txt'
export DOWNLOAD_PATH="/home/pdutta/Data/Cancer_wiseGDC/Data/"$CANCER_TYPE"/BAM/"
export VCF_PATH="/home/pdutta/Data/Cancer_wiseGDC/Data/"$CANCER_TYPE"/VCF/"

chmod 600 $TOKEN_PATH

cd ..
python GDC_client_BAM2VCF_hybrid.py $TOKEN_PATH $REF_PATH $MANIFEST_PATH $DOWNLOAD_PATH $VCF_PATH