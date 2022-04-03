cd /projects/e31039/MK/STEPS/
module load sratoolkit

###to download SRA file
prefetch SRR123456 SRR768910 # could be one or multiple SRA file


## to download SRA files and convert them to fastq format
##use --split-files  for paired end
fastq-dump --outdir ./fastq/ --split-files SRR123456


