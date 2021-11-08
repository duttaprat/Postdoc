import os, tarfile
import pyBigWig
import pybedtools
from pybedtools import BedTool
import zipfile, gzip
from gtfparse import read_gtf
import pandas as pd
pd.set_option('display.max_columns', None)


reference_genome_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/genome_by_sequence/"
reference_genome_name = "chr1.fa"
data_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/gencode_annotation/"
file_name="gencodeV38.bb"
gencode_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/gencode_annotation/"
gtf_annotation_filename= "gencode.v38.annotation.gtf"
gff_annotation_filename= "gencode.v38.annotation.gff3"
file_path= data_path+file_name
gencode_gtf_file_path= gencode_path+ gtf_annotation_filename
gencode_gff_file_path= gencode_path+ gff_annotation_filename
reference_genome_file_path = reference_genome_path + reference_genome_name
data_save_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/region_wise_sequence"


fasta = pybedtools.example_filename(reference_genome_file_path)
df1 = read_gtf(gencode_gtf_file_path)
print(df1.columns)
print(df1.shape)

df_chr1 = df1[df1['seqname']=='chr1']


df_transcript = pd.DataFrame(columns= ['transcript_id','gene_id','transcript_type','transcript_name','sequence'])
df_exon = pd.DataFrame(columns= ['exon_id','exon_number','transcript_id','gene_id','transcript_type','gene_type','sequence'])
df_gene = pd.DataFrame(columns= ['gene_id','gene_type','gene_name','sequence'])
transcript_index=exon_index = gene_index= 0
for index, row in df_chr1.iterrows():
    #print ("###", index)
    if (row['feature'] =="gene"):
        df_gene.loc[gene_index, ['gene_id']]= row['gene_id']
        df_gene.loc[gene_index, ['gene_type']]= row['gene_type']
        df_gene.loc[gene_index, ['gene_name']]= row['gene_name']
        query= row['seqname']+':'+str(row['start'])+'-'+str(row['end'])
        print (index, query)
        df_gene.loc[gene_index, ['sequence']] = BedTool.seq(query, fasta)
        gene_index= gene_index+1
    elif (row['feature'] =="transcript"):
        df_transcript.loc[transcript_index, ['transcript_id']]= row['transcript_id']
        df_transcript.loc[transcript_index, ['gene_id']]= row['gene_id']
        df_transcript.loc[transcript_index, ['transcript_type']]= row['transcript_type']
        df_transcript.loc[transcript_index, ['transcript_name']]= row['transcript_name']
        query= row['seqname']+':'+str(row['start'])+'-'+str(row['end'])
        print (index, query)
        df_transcript.loc[transcript_index, ['sequence']] = BedTool.seq(query, fasta)
        transcript_index= transcript_index+1
    elif (row['feature'] =="exon"):
        df_exon.loc[exon_index, ['exon_id']]= row['exon_id']
        df_exon.loc[exon_index, ['exon_number']]= row['exon_number']
        df_exon.loc[exon_index, ['transcript_id']]= row['transcript_id']
        df_exon.loc[exon_index, ['gene_id']]= row['gene_id']
        df_exon.loc[exon_index, ['transcript_type']]= row['transcript_type']
        df_exon.loc[exon_index, ['gene_type']]= row['gene_type']
        query= row['seqname']+':'+str(row['start'])+'-'+str(row['end'])
        print (index, query)
        df_exon.loc[exon_index, ['sequence']] = BedTool.seq(query, fasta)
        exon_index= exon_index+1


df_exon.to_csv(data_save_path +"/"+ row['seqname']+"_exon_sequence.csv", index=False)
df_transcript.to_csv(data_save_path +"/"+ row['seqname']+"_transcript_sequence.csv", index=False)
df_gene.to_csv(data_save_path +"/"+ row['seqname']+"_gene_sequence.csv", index=False)








