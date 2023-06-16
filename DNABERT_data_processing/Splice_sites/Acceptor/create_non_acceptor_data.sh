chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

for i in ${chr[@]};
do
    python create_non_acceptor_data.py $i> /home/pdutta/Data/Non_acceptor_splice_sites/Data/$i.txt
    echo "coordinates for" $i "created!!!"
done