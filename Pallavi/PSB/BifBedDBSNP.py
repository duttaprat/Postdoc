import pyBigWig
import pandas as pd


output_path= "/home/pdutta/Data/DBSNP/Chromosome_wise/"
bw = pyBigWig.open("/home/pdutta/Data/DBSNP/dbSnp155.bb")



chromosomes = bw.chroms().keys()
print(chromosomes)


for chrom in chromosomes:
    # Extract intervals for the chromosome
    if(len(chrom.split('_'))==1):
        print(chrom)
        entries = bw.entries(chrom, 0, bw.chroms()[chrom])

        if entries:  # Check if there are any entries
            # Create a list to hold the data for this chromosome
            regions = []
            for entry in entries:
                start, end, value = entry
                fields = value.split('\t')
                regions.append([chrom, start, end]+ fields)

            # Convert the regions to a DataFrame
            df = pd.DataFrame(regions, columns=['chromosome', 'start', 'end', 'name','ref','altCount','alts','shiftBases','freqSourceCount','minorAlleleFreq','majorAllele','minorAllele',	'maxFuncImpact','class','ucscNotes','_dataOffset','_dataLen'])

        # Save the dataframe to a CSV file with the chromosome name
        df.to_csv(output_path+f"{chrom}_data.csv", index=False)

# Close the BigBed file
bw.close()


# In[ ]:




