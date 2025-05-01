###########################################################################################
#      DIFFERENTIAL EXPRESSION ANALYSIS DIPLOTAXODON LIMNOTHRISSA VS BIGEYE COMPLEX       #
###########################################################################################

# Script for differential expression analysis using DESeq2 on RNA-seq data
# Last updated: 16-11-2023

# CONFIGURATION - Edit paths to match own directory structure
count_matrix_fn <- "../data/06_rnaseq_count-matrix_htseq-count.txt" 
metadata_fn <- "../data/06_rnaseq_DESeq2_metadata.txt"
ann_fn <- "../data/06_fAstCal1/genomic_biotypes.txt"
results_dir <- "path_to_DEAresults_dir"

###### Install required packages #####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install necessary packages if not already present
required_packages <- c("DESeq2", "apeglm", "RColorBrewer", "pheatmap", "ggplot2")
for (i in required_packages) {
  if (!requireNamespace(i, quietly = TRUE))
    BiocManager::install(i)}

# Load required libraries
library(DESeq2)
library(apeglm)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

packageVersion("DESeq2")
packageVersion("apeglm")

vignette("DESeq2")            # USEFUL INFO
#browseVignettes("biomaRt")    # http://127.0.0.1:31332/library/biomaRt/doc/accessing_ensembl.html


#========================================================================#
#                         DATA IMPORT                                    #
#========================================================================#

# Load data from htseq-count output and metadata (species info per sample)
data <- read.table(count_matrix_fn, header=T, row.names=1)         # gene count matrix
meta <- read.table(metadata_fn, header=T, row.names=1)             # metadata - species info

length(rownames(data))    # 26211 genes in the matrix table

# Load table containing biotypes of genes in the fAstCal1.2 genomic annotation used for the alignment of rnaseq reads.
ann <- read.csv(ann_fn, sep='\t', header=TRUE)                   # 31331 genes
ann_protein_coding <- subset(ann, BioType == 'protein_coding')   # 26070 protein-coding genes


#========================================================================#
#                     DATA FILTERING AND QC                              #
#========================================================================#

# Filter count matrix to retain only protein-coding genes
data_prot_coding <- data[which(rownames(data) %in% ann_protein_coding$GeneID),]
nrow(data_prot_coding)

# out of 26211 in the count matrix, 26070 are protein-coding.

## Create DESeq2 object
deseqObj <- DESeqDataSetFromMatrix(countData = data_prot_coding, 
                                   colData = meta, 
                                   design = ~ group)   # group = limno (limno complex), bigeye (bigeye complex)

## QC of samples prior to DE analysis
# PCA
vsd <- vst(deseqObj, blind=FALSE)    # Transform data using variance stabilizing transformation
# Blind=FALSE is set so transformation is not completely blind to sample info
head(assay(vsd), 3)

# PCA by species
pca <- plotPCA(vsd, intgroup="species", returnData = F)
pca + geom_label(aes(label = rownames(meta)))
# Option to save plot:
# ggsave(filename = "results_dir/PC1_vs_PC2_by_species.png", plot = pca, device = "png")

# PCA by group
pca <- plotPCA(vsd, intgroup="group", returnData = F)
pca + geom_label(aes(label = rownames(meta)))
# Option to save plot:
# ggsave(filename = "results_dir/PC1_vs_PC2_by_group.png", plot = pca, device = "png")


# D. holochromis sample D34G04 appears as an outlier --> exclude from analysis. The rest of 'limnothrissa' samples are retained.

# Exclude holochromis sample:
data <- data[, -which(names(data) == "D34G04")]                  # exclude column with sample from count matrix
meta <- meta[row.names(meta) != "D34G04", , drop = FALSE]        # and from metadata file

# checking if column names in data correspond to metadata samples
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# Re-create the filtered dataset and DESeq object
data_prot_coding <- data[which(rownames(data) %in% ann_protein_coding$GeneID),]
deseqObj <- DESeqDataSetFromMatrix(countData = data_prot_coding, 
                                   colData = meta, 
                                   design = ~ group)

## Filter object to remove low-expressed genes (<5 counts in less than 3 samples)
smallestGroupSize <- 3
keep <- rowSums(counts(deseqObj) >= 5) >= smallestGroupSize
deseqObj <- deseqObj[keep,] 

deseqObj

# out of 26211 in the count matrix, 26070 are protein-coding. After filtering for low read counts, 20299 genes.

#========================================================================#
#                 DIFFERENTIAL EXPRESSION ANALYSIS                       #
#========================================================================#

deseqObj <- DESeq(deseqObj)

# Information about the model fit
sizeFactors(deseqObj)   # size factor per sample
colSums(counts(deseqObj))  # total number of read counts per sample
colSums(counts(deseqObj, normalized=T))  # normalized counts per sample
plotDispEsts(deseqObj)

# Build results table
contrast <- c("group", "limno", "bigeye")   # specifies the contrast tested
DEresults <- results(deseqObj, contrast=contrast, alpha = 0.05)
head(DEresults)
mcols(DEresults)$description

# Shrink log2 fold changes (recommended for visualization and ranking)
?lfcShrink

# Get coefficient name for lfcShrink function
resultsNames(deseqObj)    # This is "group_limno_vs_bigeye"

# Shrink log2 fold changes 
DEresults_shrunk <- lfcShrink(deseqObj, coef=2)    # type default "apeglm"
head(DEresults_shrunk)

# Summary of differential expression analysis results
summary(DEresults)        # 20299  genes with non-zero value (included in analysis)
summary(DEresults_shrunk)

sum(DEresults$padj < 0.05, na.rm=TRUE)          # 4648 genes DE (before logFC shrinkage)
sum(DEresults_shrunk$padj < 0.05, na.rm=TRUE)   # 4648 genes DE at Adjpval<0.05
sum(DEresults_shrunk$padj < 0.01, na.rm=TRUE)   # 2743 genes DE at Adjpval<0.01

#========================================================================#
#                 DEA RESULTS EXPORT                                     #
#========================================================================#

# Add gene IDs as a column for easier downstream processing
DEresults_shrunk$entrezgene_accession <- rownames(DEresults_shrunk)

# Export results table
write.table(DEresults_shrunk, 
            file = "results_dir/DESeq_output_limno-vs-bigeye.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
