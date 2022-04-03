#!/bin/bash
#SBATCH -A e31039 # Allocation
#SBATCH -p short # Queue
#SBATCH -t 04:00:00 # Walltime/duration of the job
#SBATCH -N 1 # Number of Nodes
#SBATCH --mem=36G # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --ntasks-per-node=6 # Number of Cores (Processors)
#SBATCH --job-name="RNAseqAssignment" # Name of job
# add a project directory to your PATH (if needed)
export PATH=$PATH:/projects/e31039/tools

# load modules you need to use
module load python/anaconda
module load fastqc
module load samtools/1.6
module load python/anaconda3.6
module load deeptools/3.1.1
module load R/3.6.2
module load STAR/2.6.0
module load rsem/1.3.0


### GEO GSE97831

# Set your working directory
cd $SLURM_SUBMIT_DIR
DATAfq="/projects/e31039/Assignments/RNAseq"
RSEMindx="/projects/e31039/data/genomes_index"

#Perform QC using FastQC on the fastq files
mkdir FASTQC
fastqc $DATAfq/ -o ./FASTQC

#Not trimming files because FASTQC did not show significant overrepresented sequences
#Alignment and expression estimates  using the STAR-RSEM pipeline. Fastq files are mate paires from paired end sequencing. For that reason using --paired-end option.
#Preparation of Genome indexes and reference for RSEM (Not run)

#STAR --runThreadN 10 --runMode genomeGenerate --genomeDir $RSEMindx/RsemStarHg38/ --genomeFastaFiles $RSEMindx/RsemStarHg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile $RSEMindx/RsemStarHg38/Homo_sapiens.GRCh38.86.chr.gtf --sjdbOverhang 99
#rsem-prepare-reference --gtf $RSEMindx/RsemStarHg38/Homo_sapiens.GRCh38.86.chr.gtf --star --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p 10 $RSEMindx/RsemStarHg38/ $RSEMindx/RsemStarHg38/hg38StRsm

rsem-calculate-expression --star --no-bam-output --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p 10 --paired-end $DATAfq/SRR5448863_1.fastq $DATAfq/SRR5448863_2.fastq $RSEMindx/RsemStarHg38/hg38StRsm Control1

rsem-calculate-expression --star --no-bam-output --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p 10 --paired-end $DATAfq/SRR5448864_1.fastq $DATAfq/SRR5448864_2.fastq $RSEMindx/RsemStarHg38/hg38StRsm Control2

rsem-calculate-expression --star --no-bam-output --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p 10 --paired-end $DATAfq/SRR5448865_1.fastq $DATAfq/SRR5448865_2.fastq $RSEMindx/RsemStarHg38/hg38StRsm BMI1KD1

rsem-calculate-expression --star --no-bam-output --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p 10 --paired-end $DATAfq/SRR5448866_1.fastq $DATAfq/SRR5448866_2.fastq $RSEMindx/RsemStarHg38/hg38StRsm BMI1KD2

rsem-calculate-expression --star --no-bam-output --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p 10 --paired-end $DATAfq/SRR5448867_1.fastq $DATAfq/SRR5448867_2.fastq $RSEMindx/RsemStarHg38/hg38StRsm RNFKD1

rsem-calculate-expression --star --no-bam-output --star-path /software/STAR/2.6.0/bin/Linux_x86_64_static -p 20 --paired-end $DATAfq/SRR5448868_1.fastq $DATAfq/SRR5448868_2.fastq $RSEMindx/RsemStarHg38/hg38StRsm RNFKD2

#Generating data matrix of expected counts
rsem-generate-data-matrix Control1.genes.results Control2.genes.results BMI1KD1.genes.results BMI1KD2.genes.results RNFKD1.genes.results RNFKD2.genes.results > output.genes.counts.matrix.CTRLBM1KD.RNFKDtxt

rsem-generate-data-matrix Control1.isoforms.results Control2.isoforms.results BMI1KD1.isoforms.results BMI1KD2.isoforms.results RNFKD2.isoforms.results RNFKD2.isoforms.results > output.isoforms.counts.matrix.CTRLBM1KD.RNFKD.txt


rsem-run-ebseq gene.counts.matrix.CTRLBM1KD 2,2,2 CTRLBM1KDRNF.DEG
rsem-run-ebseq --ngvector  $RsemStarHg38/hg38StRsm.ngvec.ngvec output.isoforms.counts.matrix.CTRLBM1KD.RNFKD.txt 2,2,2 CTRLBM1KDRNF.DEI



