import os, tarfile
import pyBigWig
import pybedtools
import zipfile, gzip
from gtfparse import read_gtf
import pandas as pd
pd.set_option('display.max_columns', None)
from Bio.Seq import Seq


reference_genome_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/genome_by_sequence/"
data_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/gencode_annotation/"
file_name="gencodeV38.bb"
gencode_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/gencode_annotation/"
gtf_annotation_filename= "gencode.v38.annotation.gtf"
gff_annotation_filename= "gencode.v38.annotation.gff3"
file_path= data_path+file_name
gencode_gtf_file_path= gencode_path+ gtf_annotation_filename
gencode_gff_file_path= gencode_path+ gff_annotation_filename
#reference_genome_file_path = reference_genome_path + reference_genome_name
data_save_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/region_wise_sequence/"


chr_files= os.listdir(reference_genome_path)
df = read_gtf(gencode_gtf_file_path)


print(df.columns)
print(df.shape)


def acceptor_sequence(entry, df_acceptor_splice, acceptor_index):
    print ("@@@@", acceptor_index)
    df_acceptor_splice.loc[acceptor_index, ['transcript_id']]= entry['transcript_id']
    df_acceptor_splice.loc[acceptor_index, ['exon_id']] = entry['exon_id']
    df_acceptor_splice.loc[acceptor_index, ['transcript_type']] = entry['transcript_type']
    df_acceptor_splice.loc[acceptor_index, ['strand']]= entry['strand']
    if (entry['strand']=='+'):
        df_acceptor_splice.loc[acceptor_index, ['coordinates']]= int(entry['start'])-1
        acceptor_query= entry['seqname']+':'+str(int(entry['start'])+acceptor_coordinates[0])+'-'+str(int(entry['start'])+acceptor_coordinates[1])
        print ("##$$$", acceptor_query)
        df_acceptor_splice.loc[acceptor_index, ['sequence']] = pybedtools.BedTool.seq(acceptor_query, fasta)
        print ("###", acceptor_index)
    elif (entry['strand']=='-'):
        df_acceptor_splice.loc[acceptor_index, ['coordinates']]= int(entry['end'])+1
        acceptor_query= entry['seqname']+':'+str(int(entry['end'])+(acceptor_coordinates[1]*(-1)))+'-'+str(int(entry['end'])+(acceptor_coordinates[0]*(-1)))
        temp_seq = pybedtools.BedTool.seq(acceptor_query, fasta)
        rev_seq = str(Seq(temp_seq).reverse_complement())
        df_acceptor_splice.loc[acceptor_index, ['sequence']] = rev_seq


# In[45]:


def donor_sequence(entry, df_donor_splice, donor_index):
    df_donor_splice.loc[donor_index, ['transcript_id']]= entry['transcript_id']
    df_donor_splice.loc[donor_index, ['exon_id']] = entry['exon_id']
    df_donor_splice.loc[donor_index, ['transcript_type']] = entry['transcript_type']
    df_donor_splice.loc[donor_index, ['strand']]= entry['strand']
    if (entry['strand']=='+'):
        df_donor_splice.loc[donor_index, ['coordinates']] = int(entry['end'])+1
        donor_query= entry['seqname']+':'+str(int(entry['end'])+donor_coordinates[0])+'-'+str(int(entry['end'])+donor_coordinates[1])
        print ("##$$$", donor_query, donor_index)
        df_donor_splice.loc[donor_index, ['sequence']] = pybedtools.BedTool.seq(donor_query, fasta)
    elif (entry['strand']=='-'):
        df_donor_splice.loc[donor_index, ['coordinates']] = int(entry['start'])-1
        donor_query= entry['seqname']+':'+str(int(entry['start'])+(donor_coordinates[1]*(-1)))+'-'+str(int(entry['start'])+(donor_coordinates[0]*(-1)))
        print ("##$$$", donor_query, donor_index)
        temp_seq = pybedtools.BedTool.seq(donor_query, fasta)
        rev_seq = str(Seq(temp_seq).reverse_complement())
        df_donor_splice.loc[donor_index, ['sequence']] = rev_seq


# In[46]:


def promoter_sequence(entry, df_promoter, promoter_index):
    df_promoter.loc[promoter_index, ['gene_id']] = entry['gene_id']
    df_promoter.loc[promoter_index, ['transcript_id']]= entry['transcript_id']
    df_promoter.loc[promoter_index, ['transcript_type']]= entry['transcript_type']
    df_promoter.loc[promoter_index, ['strand']] = entry['strand']
    df_promoter.loc[promoter_index, ['exon_id']]=  entry['exon_id']
    if (entry['strand']=='+'):
        df_promoter.loc[promoter_index,['TSS']]= entry['start']
        promoter_query= entry['seqname']+':'+str(int(entry['start'])+promoter_coordinates[0])+'-'+str(int(entry['start'])+promoter_coordinates[1])
        df_promoter.loc[promoter_index, ['sequence']] = pybedtools.BedTool.seq(promoter_query, fasta)
    elif(entry['strand']=='-'):
        df_promoter.loc[promoter_index,['TSS']]= entry['end']
        promoter_query= entry['seqname']+':'+str(int(entry['end'])+(promoter_coordinates[1]*(-1)))+'-'+str(int(entry['end'])+(promoter_coordinates[0]*(-1)))
        temp_seq = pybedtools.BedTool.seq(promoter_query, fasta)
        #print ("Original", temp_seq)
        rev_seq = str(Seq(temp_seq).reverse_complement())
        #print (rev_seq)
        df_promoter.loc[promoter_index, ['sequence']] = rev_seq


# In[47]:


def core_promoter_sequence(entry, df_core_promoter, core_promoter_index):
    df_core_promoter.loc[core_promoter_index, ['gene_id']] = entry['gene_id']
    df_core_promoter.loc[core_promoter_index, ['transcript_id']] = entry['transcript_id']
    df_core_promoter.loc[core_promoter_index, ['transcript_type']] = entry['transcript_type']
    df_core_promoter.loc[core_promoter_index, ['strand']] = entry['strand']
    df_core_promoter.loc[core_promoter_index, ['exon_id']]  =  entry['exon_id']
    if (entry['strand']=='+'):
        df_core_promoter.loc[core_promoter_index,['TSS']]= entry['start']
        core_promoter_query= entry['seqname']+':'+str(int(entry['start'])+core_promoter_coordinates[0])+'-'+str(int(entry['start'])+core_promoter_coordinates[1])
        df_core_promoter.loc[core_promoter_index, ['sequence']] = pybedtools.BedTool.seq(core_promoter_query, fasta)
    elif(entry['strand']=='-'):
        df_core_promoter.loc[core_promoter_index,['TSS']]= entry['end']
        core_promoter_query= entry['seqname']+':'+str(int(entry['end'])+(core_promoter_coordinates[1]*(-1)))+'-'+str(int(entry['end'])+(core_promoter_coordinates[0]*(-1)))
        temp_seq = pybedtools.BedTool.seq(core_promoter_query, fasta)
        rev_seq = str(Seq(temp_seq).reverse_complement())
        df_core_promoter.loc[core_promoter_index, ['sequence']] = rev_seq





global df_promoter, df_core_promoter, df_donor_splice, df_acceptor_splice


for file in chr_files:
    if file.endswith('.fa'):
        reference_genome_file_path = reference_genome_path + file
        chr_name = file.split('.')[0]
        print (chr_name)
        df_chr= df[df['seqname']==chr_name]
        df_new = df_chr[df_chr['feature']=='exon'].reset_index(drop=True)
        print (df_new.shape)
        #print (reference_genome_file_path)
        fasta = pybedtools.example_filename(reference_genome_file_path)
        
        
        folder_name = data_save_path + chr_name   
        print (folder_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        
        df_promoter = pd.DataFrame(columns= ['gene_id','transcript_id','strand','exon_id','transcript_type','TSS', 'sequence'])
        df_core_promoter = pd.DataFrame(columns= ['gene_id','transcript_id','strand','exon_id','transcript_type','TSS','sequence'])
        df_donor_splice = pd.DataFrame(columns= ['transcript_id','exon_id','transcript_type','strand','coordinates','sequence'])
        df_acceptor_splice= pd.DataFrame(columns= ['transcript_id','exon_id', 'transcript_type','strand','coordinates', 'sequence'])
        promoter_index = core_promoter_index = donor_index = acceptor_index = 0
        
        
        
        
        df_by_transcript = df_new.head(20).groupby("transcript_id")
        for transcript, entries in df_by_transcript:
            print ("\n", transcript,len( entries) )
            print (entries.iloc[0]['strand'])
            if (entries.iloc[0]['strand']=='+'):
                for index, entry in entries.iterrows():
                    if (int(entry['exon_number'])==1):
                        print("TSS", entry['start'])
                        promoter_sequence(entry, df_promoter, promoter_index)
                        promoter_index= promoter_index+1
                        core_promoter_sequence(entry, df_core_promoter, core_promoter_index)
                        core_promoter_index = core_promoter_index +1
                        donor_sequence(entry, df_donor_splice, donor_index)
                        donor_index =donor_index+ 1
                        print("first donor_splice_site", entry['end'])
                    elif (int(entry['exon_number'])==len(entries)):
                        print ("Last acceptor_splice_site",  entry['start'])
                        acceptor_sequence(entry, df_acceptor_splice, acceptor_index)
                        acceptor_index = acceptor_index +1
                    else:
                        print ("Acceptor splice site", entry['start'])
                        acceptor_sequence(entry, df_acceptor_splice, acceptor_index)
                        acceptor_index = acceptor_index +1
                        print ("Donor splice site", entry['end'])
                        donor_sequence(entry, df_donor_splice, donor_index)
                        donor_index =donor_index+ 1
            elif (entries.iloc[0]['strand']=='-'):
                for index, entry in entries.iterrows():
                    if (int(entry['exon_number'])==1):
                        print("TSS", entry['end'])
                        print("first donor_splice_site", entry['start'])
                        promoter_sequence(entry, df_promoter, promoter_index)
                        promoter_index= promoter_index+1
                        core_promoter_sequence(entry, df_core_promoter, core_promoter_index)
                        core_promoter_index = core_promoter_index +1
                        donor_sequence(entry, df_donor_splice, donor_index)
                        donor_index =donor_index+ 1
                    elif (int(entry['exon_number'])==len(entries)):
                        print ("Last acceptor_splice_site",  entry['end'])
                        acceptor_sequence(entry, df_acceptor_splice, acceptor_index)
                        acceptor_index = acceptor_index +1
                    else:
                        print ("Acceptor splice site", entry['end'])
                        acceptor_sequence(entry, df_acceptor_splice, acceptor_index)
                        acceptor_index = acceptor_index +1
                        print ("Donor splice site", entry['start'])
                        donor_sequence(entry, df_donor_splice, donor_index)
                        donor_index =donor_index+ 1
        df_promoter.to_csv(folder_name +"/"+ entry['seqname']+"_promoter.csv", index=False)
        df_core_promoter.to_csv(folder_name +"/"+ entry['seqname']+"_core_promoter_sequence.csv", index=False)
        df_acceptor_splice.to_csv(folder_name +"/"+ entry['seqname']+"_acceptor_splice_set_sequence.csv", index=False)
        df_donor_splice.to_csv(folder_name +"/"+ entry['seqname']+"_donor_splice_set.csv", index=False) 