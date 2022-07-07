import pandas as pd
import sys
chr_name = sys.argv[1]


file_path = "/data05/shared/pdutta_data/DNABERT_2.0/pretrain/chromosome_wise_NSP_sequence/"+chr_name+"_NSP_sequence.csv"
data_save_path = "/data05/shared/pdutta_data/DNABERT_2.0/pretrain/chromosome_wise_kmer_paragraph"
df_file = pd.read_csv(file_path)


def seq2kmer(seq, k):
    seq=seq.upper()
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

grouped_by_transcript = df_file.groupby('transcript_id', sort=False)
paragraph_list = []
count = 0
print ("Number of transcript: ",len(grouped_by_transcript.size()))
for name, group in grouped_by_transcript:
    count = count+1
    if ((count%100)==0):
        print(count)
    temp_list = []
    group = group.sort_values(by=['start']).reset_index(drop=True)
    for index, row in group.iterrows():
        temp_list.append(seq2kmer(row['sequence'], 6))
    paragraph_list.append(temp_list)




file_path = data_save_path+"/"+chr_name+"_kmer_paragraph.txt"
with open(file_path, 'w') as fp:
    fp.write("[")
    for item in paragraph_list:
        fp.write("%s," % item)
    fp.write("]")






