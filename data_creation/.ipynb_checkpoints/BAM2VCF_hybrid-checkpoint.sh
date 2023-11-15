chmod +x /home/pdutta/gdc-client
export PATH=$PATH:/home/pdutta
chmod 600 /home/pdutta/Data/GDC/RAMANA/Token/gdc-user-token.2023-08-18T17_34_29.807Z.txt

python GDC_client_BAM2VCF_hybrid.py /home/pdutta/Rekha_LabWork/Collaborator_work/RNA-seq/RNAseq_CommonFiles/STAR/Genome/GRCh38.primary_assembly.genome.fa bam_manifest_file_4_2.txt 