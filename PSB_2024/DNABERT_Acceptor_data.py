import os, glob
import pandas as pd
import pysam
import numpy as np
from Bio import SeqIO



# In[2]:

chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
chr_wise_files = glob.glob('/data/private/pdutta/PSB_Data_2024/Acceptor/Chrwise/*.tsv')
reference_path  = "/data/projects/Resources/HumanReferenceGenome/"
output_path = "/data/projects/PSB/DNABERT_data/Acceptor/" 


# In[3]:


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

def get_sequences(df, genome):
    data = []
    for idx, row in df.iterrows():
        #print(row)
        if((idx %100000) ==0):
            print(row[0], idx)
        chrom = row[0]
        ref_start = row[1]
        ref_end = row[2]
        variant_start = row[11]
        variant_end = row[12]
        ref_nucleotide = row[14]
        alts = row[16].rstrip(',').split(',')
        
        # Get reference sequence
        ref_seq = str(genome[chrom].seq[ref_start:ref_end])
        #print("reference sequence" , ref_seq)

        # Adjust for 0-based indexing in python
        variant_pos_start = variant_start - ref_start
        variant_pos_end = variant_end - ref_start

        # Replace reference nucleotide with each alternative
        for alt in alts:
            if alt == "":  # Deletion
                # Fetch additional bases from genome to maintain sequence length
                extra_bases = str(genome[chrom].seq[ref_end:ref_end + len(ref_nucleotide)])
                alt_seq = ref_seq[:variant_pos_start] + ref_seq[variant_pos_end:] + extra_bases
            else:  # SNPs, Insertions, and partial Deletions
                extra_bases = str(genome[chrom].seq[ref_end:ref_end + len(ref_nucleotide) - len(alt)])
                alt_seq = ref_seq[:variant_pos_start] + alt + ref_seq[variant_pos_end:] + extra_bases
                
            # Append to the list as a dictionary
            data.append({
                'chr': row[0],
                'strand': row[4],
                'Acceptor_coordinates': str(row[1])+'-'+str(row[2]),
                'Ensemble_Transcript_ID': row[8],
                'Transcript_coordinates': str(row[6])+'-'+str(row[7]),
                'dbsnp_id': row[13],
                'varinat_start': variant_start,
                'variant_end': variant_end,
                'ref_neucleotide': ref_nucleotide,
                'alternative_neucleotide': alt,
                'reference_seq': ref_seq,
                'alt_seq': alt_seq
            })
            #print(row[13],row[8], row[0], variant_start, variant_end, ref_seq, alt_seq)
            #input()
            
    
    # Convert the list of dictionaries to a DataFrame
    new_df = pd.DataFrame(data)
    print(new_df.shape)
    new_df = new_df.drop_duplicates().reset_index()
    print(new_df.shape)
    groupby_columns = [col for col in new_df.columns if col != 'Ensemble_Transcript_ID']
    new_df = new_df.groupby(groupby_columns)['Ensemble_Transcript_ID'].apply(','.join).reset_index()
    print(new_df.shape)
    data = []
    merged_list = list(zip(new_df['reference_seq'], new_df['alt_seq']))
    merged_list = [item.upper() for tup in merged_list for item in tup]
    #print(merged_list)
    kmer_lst = list(map(seq2kmer, merged_list))
    df_kmer = pd.DataFrame(kmer_lst, columns=['Sequence'])
    
#     values = [0] * (len(df_kmer) // 2) + [1] * (len(df_kmer) // 2)

#     # If the DataFrame has an odd number of rows, add one more 0 or 1 to make the length match
#     if len(df_kmer) % 2:
#         values += [np.random.choice([0, 1])]
    df_kmer['Label'] = np.random.choice([0, 1], size=len(df_kmer))
    
    print(new_df.head())
    print("&*")
    print(df_kmer.head())
    
    return new_df, df_kmer



for chromosome in chromosomes:
    # Read the csv file into a DataFrame
    chromosome_name = 'chr'+str(chromosome)
    if (chromosome_name != 'chrM'):
        print(chromosome_name)
        genome = SeqIO.to_dict(SeqIO.parse(reference_path+chromosome_name+".fa", "fasta"))
        df = pd.read_csv('/data/private/pdutta/PSB_Data_2024/Acceptor/Chrwise/intersected_DBSNP_transcript_'+chromosome_name+'_data.tsv', header=None, sep= '\t')
        print(df.shape)
        #df= df[df[15]==2]
        #df = df.head(15)
        print(df.head())
        new_df, df_kmer = get_sequences(df, genome)
        
        out_folder_path = output_path+chromosome_name
        if not os.path.exists(out_folder_path):
            os.makedirs(out_folder_path)
        new_df.to_csv(out_folder_path + "/all_data_new.tsv", sep="\t", index= False)
        df_kmer.to_csv(out_folder_path + "/dev.tsv", sep="\t", index= False)
        print("All the files of ", chromosome_name, "are saved !!!\n ")

