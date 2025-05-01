# This script computes the integrated haplotype score (iHS) for a population ('pop') from VCF files, using the R package rehh.
# Here, this is done separately per chromosome.
# For more information, see: https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html#the-ihs-within-population-statistic

library(rehh)
library(vcfR)
library(data.table)

# using args to read input files from slurm submission script - see example in "05_rehh_xpehh_jobsubmission.slm". 
args=commandArgs(trailingOnly=TRUE)

pop <- args[1]       # Population name
outdir <- args[2]    # Output directory
chrom <- args[3]     # Chromosome
hap_fn = args[4]     # Input VCF (e.g., "${chrom}.${pop}.vcf.gz")

ihs_scan <- function(hhfile) {                      # function to calculate iHS from haplofile
    hhscan <- rehh::scan_hh(hhfile, threads = 16)
    ihh <- rehh::ihh2ihs(hhscan)
    return(ihh)
}
                              
hh <- rehh::data2haplohh(hap_file = hap_fn,              # haplotype file
                         polarize_vcf = FALSE, 
                         chr.name = chrom,
                         vcf_reader = "vcfR")

ihs <- ihs_scan(hh)

write.csv(ihs$ihs, paste(outdir,"ihs_",chrom, ".", pop, ".csv", sep = ""))
