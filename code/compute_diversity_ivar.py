# coding: utf-8 
#!/usr/bin/env python

#import the packages I need
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from collections import Counter


#import the arguments of the function
parser = argparse.ArgumentParser(description="Compute diversity scores from a vcf file computed with IVAR")
parser.add_argument("vcf_file", help = "file with the frequencies of variants per position")
parser.add_argument("downsampling", help="threshold coverage for downsampling, 90percent of the sites should be above", type=int)
parser.add_argument("frequency_threshold", help="threshold for the frequency from which we look at the variants", type = float)
parser.add_argument("length", help="length of genome", type = int)
parser.add_argument("output", help="outputfile for dataframe with the values")


args = parser.parse_args()

#threshold for downsampling
downsampling = args.downsampling
freq = args.frequency_threshold
length = args.length
output = args.output

#import vcf file as a dataframe and make sure columns are in good format
VCF_df = pd.read_csv(args.vcf_file, sep="\t")
VCF_df["POS"].astype(int)
VCF_df["ALT_DP"].astype(int)
VCF_df["TOTAL_DP"].astype(int)

min_coverage = 500


#keep only the variants that passed and that are not indels 
VCF_noindels_df = VCF_df[VCF_df.ALT.isin(["A","T","C","G"])]
VCF_noindels_pass_df = VCF_noindels_df[VCF_noindels_df.PASS]
VCF_noindels_pass_thresh_df = VCF_noindels_pass_df[VCF_noindels_pass_df.TOTAL_DP > min_coverage]

### Look only at those for which the frequency is more that the threshold "freq" ~1%
Study_df = VCF_noindels_pass_thresh_df[VCF_noindels_pass_thresh_df.ALT_FREQ >= freq]

# downsampling ~5000 
if np.percentile(Study_df["TOTAL_DP"], 90) < downsampling: 
    print("!!!!!!WARNING!!!!!! Depth of file is too low for downsampling: \n 90 percent of the loci have a coverage of "+str(np.percentile(VCF_df["TOTAL_DP"], 90)))


class Compute_diversity: 
    def several_variants(self, df_table):
        several_var = [i for i in dict(Counter(df_table.POS)) if dict(Counter(df_table.POS))[i] > 1]
        if len(several_var) > 0:
            several_variants_info = []
            for pos in several_var:
                variants =  list(df_table[df_table.POS == pos].ALT.values)
                variants_depth = list(df_table[df_table.POS == pos].ALT_DP.values)
                variants_freq = list(df_table[df_table.POS == pos].ALT_FREQ.values)
                ref_info = list(df_table[df_table.POS == pos][["REGION","POS", "REF", "REF_DP", "TOTAL_DP"]].values[0])
                ref_info.extend([variants, variants_depth, variants_freq])
                several_variants_info.append(ref_info)
            
            several_variants_df = pd.DataFrame(several_variants_info, columns= ["REGION","POS", "REF", "REF_DP", "TOTAL_DP",  "ALT", "ALT_DP", "ALT_FREQ"])
            To_keep_final = df_table[~df_table.POS.isin(several_var)][["REGION","POS", "REF", "REF_DP", "TOTAL_DP",  "ALT", "ALT_DP", "ALT_FREQ"]]
            df_table = pd.concat([To_keep_final, several_variants_df])
        #list of positions for which there are several variants.
        df_table = df_table.assign(ALT_FREQ=lambda df: [[x] if type(x) != list else x for x in df_table.ALT_FREQ])
        #sort per position
        df_table.sort_values("POS", inplace=True)
        #compite freq of alternative base
        for x in df_table.ALT_FREQ:
            x.append(1-sum(x))
        #downsampling
        df_table = df_table.assign(DOWN_SAMPLING=lambda x: [np.random.multinomial(downsampling, x, size=1)[0] for x in df_table.ALT_FREQ])
        #proba from sownsampling
        df_table = df_table.assign(DOWN_PROBA=lambda df: [[int(i)/downsampling for i in x] for x in df_table.DOWN_SAMPLING])
        #proba depending on the coverage
        df_table = df_table.assign(PROBA=[self._which_proba(x,y,z) for x,y,z in zip(df_table.TOTAL_DP, df_table.ALT_FREQ, df_table.DOWN_PROBA)])
        #entropy and NT diversity
        df_table = df_table.assign(SHANNON_ENTROPY= [self._Shannon_entropy(x) for x in df_table.PROBA])
        df_table = df_table.assign(NUCLEOTIDE_DIVERSITY= [self._Nucleotide_diversity(x, downsampling) for x in df_table.DOWN_SAMPLING])
        df_table = df_table.assign(EXPECTED_NT_DV= [self._Expected_Nucleotide_diversity(x) for x in df_table.PROBA])
        return(df_table)

   
    # if the coverage of the site is less than downsampling we don't change the ALT nb otherwise we take the downsampled one
    def _which_proba(self, depth_site,obs_proba,sampled_proba):
        if int(depth_site) > downsampling:
            proba = sampled_proba
        else:
            proba = obs_proba
        return(proba)

    def _Shannon_entropy(self, probas):
        S = 0
        for p in probas:
            if p > 0:
                S += p*np.log(p)
        return(abs(S))
    
    def _Nucleotide_diversity(self, alleles, depth):
        D = 0
        for n in alleles:
            if n > 0:
                D += n*(n-1)
        return((depth*(depth-1)-D)/(depth*(depth)-1))
    
    def _Expected_Nucleotide_diversity(self, proba):
        E_D = 0
        for p in proba:
            if p > 0:
                E_D += p**2
        return(1-E_D)
    

list_genes = set(Study_df.REGION)
#print(list_genes)
appended_data = []
for gene in list_genes:
    Study_df_gene = Study_df[Study_df.REGION == gene]
    To_keep = Study_df_gene[[ "REGION","POS", "REF", "REF_DP", "TOTAL_DP",  "ALT", "ALT_DP", "ALT_FREQ"]]
    Diversity_df = Compute_diversity().several_variants(To_keep)
    appended_data.append(Diversity_df)
    #print(appended_data)
    
Final_df = pd.concat(appended_data)


print("Nucleotide_diversity", sum(Final_df.NUCLEOTIDE_DIVERSITY)/length)
print("Expected Nucleotide_diversity", sum(Final_df.EXPECTED_NT_DV)/length)
print("Shannon_entropy",sum(Final_df.SHANNON_ENTROPY)/length )

Final_df.to_csv(output, sep="\t", index=False)
