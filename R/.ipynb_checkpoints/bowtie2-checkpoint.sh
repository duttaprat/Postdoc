GenomeIndexName=/home/pdutta/ncbi/reference/GRCh38_noalt_as/GRCh38_noalt_as
FASTAQ_1=/home/pdutta/ncbi/fastq/SRR5448863_1.fastq
FASTAQ_2=/home/pdutta/ncbi/fastq/SRR5448863_2.fastq
OUT_DIR=/home/pdutta/ncbi/bowtie_output
mkdir -p -- "$OUT_DIR"
bowtie2 --no-unal -p 20 -x $GenomeIndexName -1 "$FASTAQ_1" -2 "$FASTAQ_2" -S "$OUT_DIR"/output.sam
#bowtie2 --no-unal -p n -x index_name -1 reads_1.fastq -2 reads_2.fastq -S output.sam
