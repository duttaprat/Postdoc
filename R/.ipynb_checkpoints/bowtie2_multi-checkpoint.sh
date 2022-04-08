GenomeIndexName=/home/pdutta/ncbi/reference/GRCh38_noalt_as/GRCh38_noalt_as
OUT_DIR=/home/pdutta/ncbi/bowtie_output
mkdir -p -- "$OUT_DIR"

for i in SRR5448863 RR5448864 SRR5448866 SRR5448868
do
   echo "Welcome $i times"
   path_1="/home/pdutta/ncbi/fastq/"$i"_pass_1.fastq"
   path_2="/home/pdutta/ncbi/fastq/"$i"_pass_2.fastq"
   echo $path_1
   echo $path_2
   bowtie2 --no-unal -p 20 -x $GenomeIndexName -1 "$path_1" -2 "$path_2" -S "$OUT_DIR"_"$i"_SAM.sam 2>&1 | tee "$OUT_DIR/$i_bowtie_output.txt"
done