chmod +x /home/pdutta/gdc-client
export PATH=$PATH:/home/pdutta


export CANCER_TYPE="Lung" 
export MANIFEST_PATH="/home/pdutta/Data/Cancer_wiseGDC/Manifest_files/BronchusLung_gdc_manifest.2023-08-30_4.txt"
export TOKEN_PATH='/home/pdutta/Data/GDC/RAMANA/Token/gdc-user-token.2023-08-18T17_34_29.807Z.txt'
export DOWNLOAD_PATH="/home/pdutta/Data/Cancer_wiseGDC/Data/"$CANCER_TYPE"/BAM_NEED_TO_BE_PROCESSED/"

chmod 600 $TOKEN_PATH

cd ..
python GDC_client2BAM.py $TOKEN_PATH $MANIFEST_PATH $DOWNLOAD_PATH