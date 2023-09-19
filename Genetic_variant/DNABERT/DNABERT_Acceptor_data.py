import os, glob
import pandas as pd
import pysam
import numpy as np
import pickle



# In[2]:

cancer_type = "Brain"
intersected_data = "/data/projects/GDC_Cancer_Wise/"+cancer_type+"/Data/merged_data.pkl"
reference_path  = "/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"
dnabert_output_path = "/data/private/pdutta/DNABERT_data/Genetic_variant/"+ cancer_type 
raw_output_path = "/data/projects/GDC_Cancer_Wise/"+cancer_type+"/Data/"



def seq2kmer(seq):
    """
    Convert original sequence to kmers
    
    Arguments:
    seq -- str, original sequence.
    k -- int, kmer of length k specified.
    
    Returns:
    kmers -- str, kmers separated by space
    """
    k=6
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers


# In[7]:

def get_sequences(df, reference_fasta):
    data = []
    for idx, row in df.iterrows():
        #print(row)
        chrom = row['chr_name']
        ref_start = row['start']
        ref_end = row['end']
        variant_start = row['START_POS']
        variant_end = row['END_POS']
        ref_nucleotide = row['REF']
        alt = row["ALT"]
        
        # Adjust for 0-based indexing in python
        variant_pos_start = variant_start - ref_start
        variant_pos_end = variant_end - ref_start
        #print(ref_nucleotide , alt, variant_pos_start, variant_pos_end)
        
        
        # Get reference sequence
        ref_seq = reference_fasta.fetch(chrom, ref_start, ref_end)
        #print(ref_seq)
        
        # Handle insertion and deletion to get the correct alt sequence
        # Identify if the variant is an insertion or deletion
        if len(ref_nucleotide) < len(alt):  # Insertion
            delete_size =  len(alt) - len(ref_nucleotide)
            #print(variant_pos_start, variant_pos_end ,delete_size)
            alt_seq = ref_seq[:variant_pos_start] + alt + ref_seq[variant_pos_end:len(ref_seq) - delete_size]

        elif len(ref_nucleotide) > len(alt):  # Deletion
            insert_size = len(ref_nucleotide) - len(alt)
            #print(insert_size)
            extra_bases = reference_fasta.fetch(chrom, ref_end, ref_end + insert_size)
            #print(extra_bases)
            alt_seq = ref_seq[:variant_pos_start] + alt + ref_seq[variant_pos_end:] + extra_bases

        else:  # SNV
            alt_seq = ref_seq[:variant_pos_start] + alt + ref_seq[variant_pos_end:]
        # print(ref_seq)
        # print(alt_seq)
        # input()


        data.append({
            'chr': chrom,
            'strand': row['strand'],
            'Transcript_ID': row['transcript_id'],
            'Acceptor_start': row['start'],
            'Acceptor_end': row['end'],
            'varinat_start': variant_start,
            'variant_end': variant_end,
            'ref_neucleotide': ref_nucleotide,
            'alternative_neucleotide': alt,
            'reference_seq': ref_seq,
            'alt_seq': alt_seq
        })
            
    
    # Convert the list of dictionaries to a DataFrame
    new_df = pd.DataFrame(data)
    # print(new_df.shape)
    # print(new_df)
    new_df = new_df.drop_duplicates().reset_index()
    print(new_df.shape)
    print("&*")
    data = []
    merged_list = list(zip(new_df['reference_seq'], new_df['alt_seq']))
    merged_list = [item.upper() for tup in merged_list for item in tup]
    #print(merged_list)
    kmer_lst = list(map(seq2kmer, merged_list))
    df_kmer = pd.DataFrame(kmer_lst, columns=['Sequence'])
    df_kmer['Label'] = np.random.choice([0, 1], size=len(df_kmer))
    print(df_kmer.shape)
    return new_df, df_kmer


## Load the VCF files from pickle
with open(intersected_data, "rb") as file:
    loaded_dictionary = pickle.load(file)
print("****All the VCF files are loaded***")   
## Load the reference file
reference_fasta = pysam.FastaFile(reference_path)
print("#####The reference file is loaded######")
# In[3]:

dnabert_raw_data = {}
for key, value in loaded_dictionary.items():
    #display(value[['chr_name', 'start', 'end', 'strand', 'transcript_type', 'CHROM', 'START_POS', 'END_POS', 'ID', 'REF', 'ALT']])
    print("CODE STARTED")
    print(key, value.shape)
    #ins_del_df = value[(value['ALT'].str.len() > value['REF'].str.len()) | (value['REF'].str.len() > value['ALT'].str.len())]
    new_df, df_kmer = get_sequences(value, reference_fasta)
    dnabert_raw_data[key.split("_")[0]] = new_df
    out_folder_path = dnabert_output_path+ "/"+ key.split("_")[0]
    print(out_folder_path, "is created")
    #input()
    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)
    df_kmer.to_csv(out_folder_path + "/dev.tsv", sep="\t", index= False)
    #print("All the files of ", chromosome_name, "are saved !!!\n ")

with open(raw_output_path+"/raw_acceptor_vcf_data.pkl", "wb") as file:
    pickle.dump(dnabert_raw_data, file)




