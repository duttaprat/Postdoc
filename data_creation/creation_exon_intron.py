import os, tarfile
import pyBigWig
import pybedtools
import zipfile, gzip
from gtfparse import read_gtf
import pandas as pd
from pandas_profiling import ProfileReport
import sys



chr_name = sys.argv[1]


# In[3]:


reference_genome_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/genome_by_sequence/"
reference_genome_name = chr_name + ".fa"
data_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/gencode_annotation/"
file_name="gencodeV38.bb"
gencode_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/gencode_annotation/"
gtf_annotation_filename= "gencode.v38.annotation.gtf"
gff_annotation_filename= "gencode.v38.annotation.gff3"
file_path= data_path+file_name
gencode_gtf_file_path= gencode_path+ gtf_annotation_filename
gencode_gff_file_path= gencode_path+ gff_annotation_filename
reference_genome_file_path = reference_genome_path + reference_genome_name
data_save_path= "/mnt/data05/shared/pdutta_data/DNABERT_2.0/pretrain/chromosome_wise_NSP_sequence"




df = read_gtf(gencode_gtf_file_path)
print(df.columns)
print(df.shape)
print (df.groupby(['feature']).size())




df = df[(df['seqname']==chr_name)].reset_index(drop=True)
df_chr = df[(df['feature']=='transcript') | (df['feature']=='exon')].reset_index(drop=True)





profile = ProfileReport(df_chr, title="GTF dataframe Profiling Report", minimal=True)
profile.to_file(chr_name+"_GTF_output.html")





fasta = pybedtools.example_filename(reference_genome_file_path)





df_transcript_paragraph = pd.DataFrame(columns= ['transcript_id','transcript_type','transcript_name','feature','start','end'])
grouped_by_transcript = df_chr.groupby('transcript_id', sort=False)
print ("Number of transcript: ",len(grouped_by_transcript.size()))
transcript_count =0
transcript_index = 0
for name, group in grouped_by_transcript:
    if (len(group)>2):
        transcript_count = transcript_count +1
        if (transcript_count%50 ==0):
            print (transcript_count)
        #print (len(group))
        group = group.sort_values(by=['start']).reset_index(drop=True)
        for index, row in group[1:].iterrows():
            if (index==1):
                df_transcript_paragraph.loc[transcript_index, ['transcript_id']]= row['transcript_id']
                df_transcript_paragraph.loc[transcript_index, ['transcript_type']]= row['transcript_type']
                df_transcript_paragraph.loc[transcript_index, ['transcript_name']]= row['transcript_name']
                df_transcript_paragraph.loc[transcript_index, ['feature']]= row['feature']
                df_transcript_paragraph.loc[transcript_index, ['start']]= row['start']
                df_transcript_paragraph.loc[transcript_index, ['end']]= row['end']
                df_transcript_paragraph.loc[transcript_index+1, ['transcript_id']]= row['transcript_id']
                df_transcript_paragraph.loc[transcript_index+1, ['transcript_type']]= row['transcript_type']
                df_transcript_paragraph.loc[transcript_index+1, ['transcript_name']]= row['transcript_name']
                df_transcript_paragraph.loc[transcript_index+1, ['feature']]= 'intron'
                df_transcript_paragraph.loc[transcript_index+1, ['start']]= row['end']+1
                transcript_index = transcript_index +2
            elif ((index +1) == len(group)):
                df_transcript_paragraph.loc[transcript_index, ['transcript_id']]= row['transcript_id']
                df_transcript_paragraph.loc[transcript_index, ['transcript_type']]= row['transcript_type']
                df_transcript_paragraph.loc[transcript_index, ['transcript_name']]= row['transcript_name']
                df_transcript_paragraph.loc[transcript_index, ['feature']]= row['feature']
                df_transcript_paragraph.loc[transcript_index, ['start']]= row['start']
                df_transcript_paragraph.loc[transcript_index, ['end']]= row['end']
                df_transcript_paragraph.loc[transcript_index-1, ['end']]= row['start']-1
                transcript_index = transcript_index+1
            else:
                df_transcript_paragraph.loc[transcript_index, ['transcript_id']]= row['transcript_id']
                df_transcript_paragraph.loc[transcript_index, ['transcript_type']]= row['transcript_type']
                df_transcript_paragraph.loc[transcript_index, ['transcript_name']]= row['transcript_name']
                df_transcript_paragraph.loc[transcript_index, ['feature']]= row['feature']
                df_transcript_paragraph.loc[transcript_index, ['start']]= row['start']
                df_transcript_paragraph.loc[transcript_index, ['end']]= row['end']
                df_transcript_paragraph.loc[transcript_index+1, ['transcript_id']]= row['transcript_id']
                df_transcript_paragraph.loc[transcript_index+1, ['transcript_type']]= row['transcript_type']
                df_transcript_paragraph.loc[transcript_index+1, ['transcript_name']]= row['transcript_name']
                df_transcript_paragraph.loc[transcript_index+1, ['feature']]= 'intron'
                df_transcript_paragraph.loc[transcript_index+1, ['start']]= row['end']+1
                df_transcript_paragraph.loc[transcript_index-1, ['end']]= row['start']-1
                transcript_index = transcript_index +2

print("Number of Transcript which more than one exon: ", transcript_count)
df_transcript_paragraph['length'] = df_transcript_paragraph['end']- df_transcript_paragraph['start']




for index, row in df_transcript_paragraph.iterrows():
    query= chr_name+':'+str(row['start'])+'-'+str(row['end'])
    df_transcript_paragraph.loc[index, ['sequence']] = pybedtools.BedTool.seq(query, fasta)
    index= index+1




df_transcript_paragraph.to_csv(data_save_path +"/"+ chr_name+"_NSP_sequence.csv", index=False)




