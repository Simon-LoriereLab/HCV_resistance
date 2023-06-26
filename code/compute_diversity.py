# coding: utf-8 
#!/usr/bin/env python

#import the packages I need
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse


#import the arguments of the function
parser = argparse.ArgumentParser(description="Compute diversity scores from a vcf file")
parser.add_argument("vcf_file", help = "file with the frequencies of variants per position")
parser.add_argument("downsampling", help="threshold coverage for downsampling, 90percent of the sites should be above", type=int)
parser.add_argument("output", help="outputfile for dataframe with the values")


args = parser.parse_args()

#threshold for downsampling
downsampling = args.downsampling
output = args.output

#import vcf file as a dataframe
VCF_df = pd.read_csv(args.vcf_file, sep="\t", comment= "#", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DEPTH"])

#Reshape the dataframe by adding new columns
VCF_df["READ_DEPTH"] = [x.split(";")[0].split("=")[-1] for x in VCF_df.INFO]

if np.percentile(VCF_df["READ_DEPTH"], 90) < downsampling:
    print("!!!!!!WARNING!!!!!! Depth of file is too low for downsampling: \n 90 percent of the loci have a coverage of "+str(np.percentile(VCF_df["READ_DEPTH"], 90)))
#remove indels positions
STUDY_DF = VCF_df[VCF_df.READ_DEPTH != "INDEL"]





STUDY_DF = STUDY_DF.assign(VAR_DEPTH=lambda df: [x.split(":")[-1] for x in STUDY_DF.DEPTH])
STUDY_DF = STUDY_DF.assign(PROBA_OBS=lambda df: [[int(i)/int(y) for i in x.split(",")] for x,y in zip(STUDY_DF.VAR_DEPTH, STUDY_DF.READ_DEPTH)])




To_keep = STUDY_DF[["POS", "REF", "ALT", "READ_DEPTH", "VAR_FREQ", "PROBA_OBS"]]

To_keep.READ_DEPTH.apply(int)



To_keep = To_keep.assign(DOWN_SAMPLING=lambda x: [np.random.multinomial(downsampling, x, size=1)[0] for x in To_keep.PROBA_OBS])

To_keep = To_keep.assign(DOWN_PROBA=lambda df: [[int(i)/downsampling for i in x] for x in To_keep.DOWN_SAMPLING])


def which_proba(depth,obs_proba,sampled_proba):
    if int(depth) > downsampling:
        proba = sampled_proba
    else:
        proba = obs_proba
    return(proba)


To_keep = To_keep.assign(PROBA=[which_proba(x,y,z) for x,y,z in zip(To_keep.READ_DEPTH, To_keep.PROBA_OBS, To_keep.DOWN_PROBA)])


def Shannon_entropy(probas):
    S = 0
    for p in probas:
        if p > 0:
            S += p*np.log(p)
    return(abs(S))
    
To_keep = To_keep.assign(SHANNON_ENTROPY= [Shannon_entropy(x) for x in To_keep.PROBA])



def Nucleotide_diversity(alleles, depth):
    D = 0
    for n in alleles:
        if n > 0:
            D += n*(n-1)
    return((depth*(depth-1)-D)/(depth*(depth)-1))
    
To_keep = To_keep.assign(NUCLEOTIDE_DIVERSITY= [Nucleotide_diversity(x, downsampling) for x in To_keep.DOWN_SAMPLING])

def Expected_Nucleotide_diversity(proba):
    E_D = 0
    for p in proba:
        if p > 0:
            E_D += p**2
    return(1-E_D)
    
To_keep = To_keep.assign(EXPECTED_NT_DV= [Expected_Nucleotide_diversity(x) for x in To_keep.PROBA])


print("Nucleotide_diversity", sum(To_keep.NUCLEOTIDE_DIVERSITY)/len(To_keep))
print("Expected Nucleotide_diversity", sum(To_keep.EXPECTED_NT_DV)/len(To_keep))
print("Shannon_entropy",sum(To_keep.SHANNON_ENTROPY)/len(To_keep) )

To_keep.to_csv(output, sep="\t", index=False)