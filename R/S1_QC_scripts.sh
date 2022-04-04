mkdir FASTQC
TRIM="/projects/e31039/softwares/Trimmomatic-0.38"
DATA="/projects/e31039/data/rnaseq"
CleanData="/projects/e31039/MK/TrimoMat/"
module load fastq
fastqc ../data/rnaseq/HBR_Rep3_chr22.R1.fastq.gz -o ./FASTQC

### Could use TrimGalore
 ## default is --phred33, --fastqc OFF -a automatic detection  --length 20
mkdir TrimGalr
module load TrimGalore
trim_galore --phred33 --paired $DATA/UHR_Rep2_chr22.R1.fastq.gz $DATA/UHR_Rep2_chr22.R2.fastq.gz --length 10 -o TrimGalr


#or use Trimmomatic
## If no quality score is specified, phred-64 is the defaul
#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clipthreshold>:<simple clip threshold><minAdapterLength>:<keepBothReads>

modulel load java 
mkdir TrimoMat
java -jar $TRIM/trimmomatic-0.38.jar PE -threads 5 -phred64 -trimlog trimM.log $DATA/HBR_Rep2_chr22.R1.fastq.gz $DATA/HBR_Rep2_chr22.R2.fastq.gz $CleanData/clean_HBR_Rep2_R1_chr22.R2.fastq.gz $CleanData/clean_HBR_Rep2_R1_unpaired.fq.gz $CleanData/clean_HBR_Rep2_chr22.R2.fastq.gz $CleanData/clean_HBR_Rep2_R2.unpaired.fq.gz ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36



