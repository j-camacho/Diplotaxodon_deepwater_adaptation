""""
This script is used to generate an empirical distribution of extreme XP-EHH scores by randomly sampling windows across the genome. 
This was used to compare XP-EHH scores in candidate regions against this background distribution.

It takes files with XP-EHH scores calculated for all chromosomes as input (see "05_rehh_xpehh_slurm.R"). Steps:
  (1) Reading of XP-EHH scores per chromosome and merging of files
  (2) Random sampling of n windows (--niter) of specified size (--winsize)
  (3) For each window, the most extreme XP-EHH score (maximum positive or negative) is extracted

The output is a text file containing a list of maximum XP-EHH scores from each random window.

Usage: 
python 05_generate_xpehh_empirical_distribution.py --winsize <window_size_in_bp> --winname <window_size_in_kb> --niter <iterations> --date <date_string>

"""

import re
import os
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-ws", "--winsize", type=int)
parser.add_argument("-wn", "--winname", type=int)
parser.add_argument("-i", "--niter", type=int)
parser.add_argument("-d", "--date", type=str)

args = parser.parse_args()

##### ARGUMENTS #####
window_size=args.winsize     # --winsize: window size in base pairs
window_size_n=args.winname   # --winname: window size in kb for file naming
niterations=args.niter       # --niter: number of random windows to sample
date=args.date               # --date: date string for output file naming

##### PATH TO FILES #####
xpehh_per_chrom="path/to/rehh_xpehh/output/"      # Directory containing rehh XP-EHH results per chromosome (can be generated using "05_rehh_xpehh_slurm.R")
outdir="path/to/output/directory/"                # Output file location

##### PRELIMINARY DATA MANIPULATION #####

# function to read input files from rehh directory in alphanumeric order
def alphanum_order(string):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string)]
    
# read in files in dir with XPEHH results per chromosome (one file per chrom)
files = sorted([fname for fname in os.listdir(xpehh_per_chrom) if fname.startswith("xpehh_chr")],key=alphanum_order)   # Adapt to match file names

# make list of dataframes to concatenate
df_list = []    

for i in range(len(files)):      # append datasets into the list
    temp_df = pd.read_csv(xpehh_per_chrom + files[i])
    df_list.append(temp_df)
    
# dataframe containing XPEHH results for all chromosomes
allchrom_df = pd.concat(df_list)
allchrom_df.rename(columns={'Unnamed: 0': 'chr_pos'}, inplace=True)   # change column name
allchrom_df['POSITION'] = allchrom_df['POSITION'].astype(int)         # change data type to int

xpehh_df = allchrom_df.set_index(['CHR'])

# Need chromosome length of all chromosomes
# Make table with this info
chrom_list = sorted(list(set(xpehh_df.index)),key=alphanum_order)   # list of chromosomes in df
chrom_lenght_list = []                                              # list of chromosome length, sorted by chromosome

for chromname in chrom_list:
    chrom_lenght_list.append(len(xpehh_df.loc[chromname]))
    
chrom_length_df = pd.DataFrame(list(zip(chrom_list,chrom_lenght_list)),   # dataframe
                               columns = ['chrom','length']).set_index('chrom')

xpehh_df = xpehh_df.reset_index().set_index(['CHR','POSITION'])


##### DISTRIBUTION OF EXTREME XPEHH VALUES IN RANDOM WINDOWS #####

empirical_scores=[]

for i in range(niterations):
    chrom = np.random.choice(chrom_list, replace=True, p=np.array(chrom_lenght_list)/np.sum(np.array(chrom_lenght_list)))
    start = np.random.choice(xpehh_df.loc[chrom].index.get_level_values('POSITION'))
    snps = xpehh_df.loc[chrom].loc[start:start+window_size]
    max_score_positive = max(snps['XPEHH_dlimno_dbigeye'])      # This should match the column containing XP-EHH scores in the input files
    max_score_negative = min(snps['XPEHH_dlimno_dbigeye'])
   
    if max_score_positive > max_score_negative:
        empirical_scores.append(max_score_positive)
    else:
        empirical_scores.append(max_score_negative)

# Save data
wf = open(f"{outdir}_{window_size_n}kb-win_{niterations}iterations_{date}.txt","w")
wf.write(str(empirical_scores))
wf.close()
                               