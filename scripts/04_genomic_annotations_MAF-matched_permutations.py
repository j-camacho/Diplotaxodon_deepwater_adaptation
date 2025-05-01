"""
This script performs an overrepresentation test for genomic annotations in GWAS outlier SNPs using MAF-matched permutations. It:
  1) Groups SNPs from the complete GWAS callset into discrete MAF bins
  2) For each iteration (n=1000), samples SNPs from these bins to match the MAF distribution of GWA outlier SNPs
  3) Counts the occurrence of each specified snpEff annotation type in the sampled SNPs
  4) Stores these counts to generate an empirical null distribution
  
Input: annotated VCF file (snpEff), MAF files for all SNPs in the GWA callset and for outlier SNPs, with columns chr', 'ps', 'maf'.
Output: Distribution of annotation counts across permutations for calculating empirical p-values.

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys 
#sys.path.append("/data/antwerpen/grp/asvardal/hs_tools")     # see https://github.com/feilchenfeldt/pypopgen3
from pypopgen3.modules import vcfpandas as vp                # https://github.com/feilchenfeldt/pypopgen3
import timeit

# Parameters - Define before running
ann_fn = "path_to_snpEff_annotated_VCF.vcf.gz"

annotations = [                          # Set annotations to test
    'intergenic_region',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'intron_variant',
    'missense_variant',
    'synonymous_variant',
    '3_prime_UTR_variant',
    '5_prime_UTR_variant',
    'splice_region_variant&intron_variant']

n_iterations = 1000
test_set = "top01pct"     # Used in output file names
date = "dd-mm-yy"
binsize = 0.05

############ Read in data ###############

# Here, "MAF files" are space-separated files with three columns: 'chr' (chromosome),'ps' (position) and 'maf' (minor allele frequency) such as:
# chr1 5623344 0.0199862

# These are provided for the top 0.01% GWA outlier SNPs and for all positions in the GWA callset

def read_snp_data(filename):
    return pd.read_csv(filename, sep=' ', names=['chr', 'ps', 'maf'])

top01pctsnp = read_snp_data("path/to/001pctGWAoutlier_MAF-per-pos.txt")
allsnp = read_snp_data("path/to/allsnp_MAF-per-pos.txt")


# Bin MAF (0.05 bins)
# MAF are grouped according to their value into different categories
def get_maf_bin(maf):
    return round(maf/0.05)

# Apply maf_bin function to datasets 
allsnp['mafbin'] = allsnp['maf'].apply(get_maf_bin)
top01pctsnp['mafbin'] = top01pctsnp['maf'].apply(get_maf_bin)

# Create a groupby object
mafgroup = allsnp.groupby('mafbin')


########## Run ################
start = timeit.default_timer()
ann_per_pos = {annotation: [] for annotation in annotations}

for i in range(n_iterations):
    intergenic_region=0
    upstream_gene_variant=0
    downstream_gene_variant=0
    intron_variant=0
    missense_variant=0
    synonymous_variant=0
    three_prime_UTR_variant=0
    five_prime_UTR_variant=0
    splice_region_variant=0
    
    for mafbin, nsnps in top01pctsnp.groupby('mafbin').apply(len).items():
        df = mafgroup.get_group(mafbin)      
        target_snps = df.sample(nsnps)
        
        for chrom,pos in target_snps[['chr','ps']].itertuples(index=False):
            header, info_dic = vp.parse_vcf_header(ann_fn.format(chrom))
            a = vp.get_vcf_df(ann_fn.format(chrom), chrom=str(chrom), start=pos ,end=pos , header=header)
            if 'intergenic_region' in a.iloc[0]['INFO'].split('ANN')[1]:
                intergenic_region += 1
            if 'upstream_gene_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                upstream_gene_variant += 1
            if 'downstream_gene_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                downstream_gene_variant += 1
            if 'intron_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                intron_variant += 1    
            if 'missense_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                missense_variant += 1   
            if 'synonymous_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                synonymous_variant += 1
            if '3_prime_UTR_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                three_prime_UTR_variant += 1
            if '5_prime_UTR_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                five_prime_UTR_variant += 1
            if 'splice_region_variant&intron_variant' in a.iloc[0]['INFO'].split('ANN')[1]:
                splice_region_variant += 1
    
    ann_per_pos['intergenic_region'].append(intergenic_region)
    ann_per_pos['upstream_gene_variant'].append(upstream_gene_variant)
    ann_per_pos['downstream_gene_variant'].append(downstream_gene_variant)
    ann_per_pos['intron_variant'].append(intron_variant)
    ann_per_pos['missense_variant'].append(missense_variant)
    ann_per_pos['synonymous_variant'].append(synonymous_variant)
    ann_per_pos['3_prime_UTR_variant'].append(three_prime_UTR_variant)
    ann_per_pos['5_prime_UTR_variant'].append(five_prime_UTR_variant)
    ann_per_pos['splice_region_variant&intron_variant'].append(splice_region_variant)

    
wf = open(f"path_to_outdir/eff-per-pos_AF-matched_{test_set}SNPs_{n_iterations}iterations_{date}.txt","w")
wf.write(str(ann_per_pos))
wf.close()

stop = timeit.default_timer()
print('Time: ', stop - start) 