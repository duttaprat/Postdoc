chr=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

for i in ${chr[@]};
do
    python create_non_donor_coordinates.py $i> /home/pdutta/Data/Splicesites_data/Non_donor_splice_sites/Coordinates/$i.txt
    echo "Data for" $i "created!!!"
done