# This script computes the cross-population extended haplotype homozygosity statistic (XP-EHH) for two given populations ('pop1', 'pop2') using the R package rehh.
# Here, VCF files are used as input, provided as separate files per population and chromosome ('hap_fn1' and 'hap_fn2' below).

# For more information see: https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html#about-the-package

# =====================================================================================================================
library(rehh)
library(vcfR)
library(data.table)

# using args to read input files from slurm submission script. As an example, see "05_rehh_xpehh_jobsubmission.slm"
args=commandArgs(trailingOnly=TRUE)

pop1 <- args[1]    # Population 1 as specified in the VCF file name (e.g., 'dlimnothrissa')
pop2 <- args[2]    # Population 2 
outdir <- args[3]  # Output directory
chrom <- args[4]   # Chromosome as specified in the VCF file name


hap_fn1 = args[5]  # VCF with population 1 samples
hap_fn2 = args[6]  # VCF population 2
                              
hhpop1 <- rehh::data2haplohh(hap_file = hap_fn1,              # haplotype file pop 1
                              polarize_vcf = FALSE, 
                              chr.name = chrom,
                              vcf_reader = "vcfR")
hhpop2 <- rehh::data2haplohh(hap_file = hap_fn2,              # haplotype file pop 2 
                              polarize_vcf = FALSE, 
                              chr.name = chrom,
                              vcf_reader = "vcfR")
    
hhscan_pop1 <- rehh::scan_hh(hhpop1, threads = 16)
hhscan_pop2 <- rehh::scan_hh(hhpop2, threads = 16)
    
xpehh <- rehh::ies2xpehh(hhscan_pop1, hhscan_pop2,
                         popname1 = pop1, popname2 = pop2,
                         include_freq = T)

# Export table. Output: one file per chromosome with XP-EHH pop1 vs pop2
write.csv(xpehh, paste(outdir,"xpehh_",chrom, ".", pop1,"_",pop2,".csv", sep = ""))
            