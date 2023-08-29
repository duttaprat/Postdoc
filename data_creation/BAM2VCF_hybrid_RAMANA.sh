chmod +x /home/campus.stonybrook.edu/pdutta/gdc-client
export PATH=$PATH:/home/campus.stonybrook.edu/pdutta

chmod 600 /data/projects/BAM_files/token/gdc-user-token.2023-08-18T17_34_29.807Z.txt

python GDC_client_BAM2VCF_hybrid.py bam_manifest_file_2_2.txt