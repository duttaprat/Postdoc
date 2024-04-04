chmod +x /home/pdutta/gdc-client

export PATH=$PATH:/home/pdutta


chmod 600 /home/pdutta/Data/Cancer_wiseGDC/GDC_token/gdc-user-token.2024-03-19T01_11_26.079Z.txt

# gdc-client download -m /home/pdutta/Data/Cancer_wiseGDC/New_data/Brain/Manifest/gdc_manifest.2024-04-04_brain.txt -t /home/pdutta/Data/Cancer_wiseGDC/GDC_token/gdc-user-token.2024-03-19T01_11_26.079Z.txt -d /home/pdutta/Data/Cancer_wiseGDC/New_data/Lung/Downloaded_files/VCF/

gdc-client download -m /home/pdutta/Data/Cancer_wiseGDC/New_data/Breast/Manifest/gdc_manifest.2024-04-04_breast.txt \
                    -t /home/pdutta/Data/Cancer_wiseGDC/GDC_token/gdc-user-token.2024-03-19T01_11_26.079Z.txt \
                    -d /home/pdutta/Data/Cancer_wiseGDC/New_data/Breast/Downloaded_files/VCF/
