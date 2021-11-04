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



df1 = read_gtf(gencode_gtf_file_path)
print(df1.columns)
print(df1.shape)


a = pybedtools.BedTool(gencode_gff_file_path)
fasta = pybedtools.example_filename(reference_genome_file_path)



a = a.sequence(fi=fasta)

print ("Crossed  the hurdle")

print(open(a.seqfn).read())

