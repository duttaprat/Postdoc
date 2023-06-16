import os, tarfile
import re, glob
import pandas as pd
import numpy as np
import zipfile, gzip
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
import plotly.graph_objects as go
from natsort import natsorted
import pybedtools
import sys

reference_genome_folder="/home/pdutta/Data/Human_Genome_Data/UCSC/genome_by_sequence/"
coordinate_path = "/home/pdutta/Data/Splicesites_data/Non_donor_splice_sites/Coordinates/"


chr_num=sys.argv[1]

def extract_substrings(sequence, pattern, padding_left, padding_right):
    n = len(sequence)
    m = len(pattern)
    substrings = []

    for i in range(n-m+1):
        if sequence[i:i+m] == pattern:
            start = max(0, i-padding_left)
            end = min(n, i+m+padding_right)
            substrings.append(sequence[start:end])

    return substrings



reference_genome_file_path = reference_genome_folder + chr_num +".fa"
bed_file_path = coordinate_path + chr_num+".bed"

fasta = pybedtools.example_filename(reference_genome_file_path)
bed = pybedtools.BedTool(bed_file_path)
bed_seq = bed.sequence(fi=fasta)
with open(bed_seq.seqfn, 'r') as f:
    for line in f:
        line = line.strip().upper()
        substrings = extract_substrings(line.strip(), "GT", 40, 38)
        for i in substrings:
            if (len(i)==80):
                print(i)