[![DOI](https://zenodo.org/badge/DOI/YOUR-DOI-HERE.svg)](https://doi.org/YOUR-DOI-HERE)


## Citation and Data Access

This repository contains the code, small data files, and notebooks associated with the study:  
**"Widespread genetic signals of visual system adaptation in deepwater cichlid fishes" [DOI:]**

### Associated Datasets

- **SNP data (VCF, per chromosome)**: [Zenodo DOI for VCF dataset](10.5281/zenodo.15576996)
- **Code and small data files (this repo)**: [Zenodo DOI for this repository](https://doi.org/...)

## Description of the data and file structure

Files and code (scripts/notebooks) are organized into:

* Data

  Raw data files associated with the study. Intermediate/output files are included when relevant for subsequent analyses.
* Scripts
* Notebooks

Files names use a numeric prefix system (e.g., 01_, 02_, etc.) to group related files by analysis. This prefix is consistent across data files, scripts, and notebooks for each analysis and the numbering sequence generally follows the order in which analyses appear in the manuscript.

Description of the files, organized by analysis:

### 1. Eye size Malawi cichlid radiation (01\_, 02\_)

#### Data
* **01_morphological_data_landmark_coordinates_Malawi_radiation.tps**

  TPS file containing coordinate data from 18 homologous landmark points on 1,376 Malawi cichlids. Eye size and body shape data were extracted from this file. See Supplementary Fig. 7 and Notebook provided ("01_extract_eye-size_data_from_TPS.ipynb") for details.

  Note: Not all specimens in this file were included in the final analysis; the subset of specimens used, along with their measurements, is provided in Supplementary Data 1.

* **01_classifier_landmark_data_morphoJ.csv**

  Classifier dataframe used to subset the TPS file for the analysis of body shape in MorphoJ (https://morphometrics.uk/MorphoJ_guide/frameset.htm?index.htm), containing:
  * im_id: photo ID
  * supplier_id: simplified specimen ID
  * spname: species name (genus_species), only shown for Diplotaxodon
  * sequenced: denotes if a specimen was sequenced (1) or not (0)
For this analysis, the TPS was filtered for 'sequenced == 1' to retain Diplotaxodon included in the GWAS for eye size (N = 53).
* **01_body_shape_diplotaxodon_GWAS-samples_PC-scores_morphoJ.txt**

  Tab-delimited file containing principal component (PC) scores from the principal component analysis (PCA) of the Procrustes coordinates for the 18 body landmark points (from the TPS file) in MorphoJ. Individual photo IDs (column 1) match the 'im_id' field in the classifier table and can be linked to sequence IDs using Supplementary Data 1.
* **02_total_pwd.biallelic_snps_pass_ancestral_as_sample.all_chrom.non_inverted.diploid.tsv**

  Table of pairwise genetic distances for 1,376 Malawi cichlids. Used to calculate genetic distance differences between Diplotaxodon individuals and to test for correlations with 1) eye size differences and 2) body shape (PC1) differences.
* **02_NJ_tree_diplotaxodon.txt**

  Neighbor-joining tree of Diplotaxodon species included in this study.

#### Scripts/Notebooks
* **01_extract_eye-size_data_from_TPS.ipynb**

  Calculation of relative eye size from landmark coordinates.
* **02_Pagels_lambda_phytools.R**

  Calculation of phylogenetic signal (Pagel's lambda) using the 'phylosig' function from the R package phytools.

### 2. Genome-wide association analysis

#### Data

* Phenotypes: relative eye size measurements for sequenced Diplotaxodon. These can be derived from Supplementary Data 1 and 3.
* Genotypes: [Zenodo DOI for VCF dataset](https://doi.org/...)

#### Workflow:

1. The VCF file was subset using bcftools to retain only Diplotaxodon samples with phenotypic data (Supplementary Data 3, column 'GWAS' == 1)
2. Binary PED files were generated using:
    `plink --vcf {VCF}`
4. Individual phenotypes (relative eye size) were added to the 6th column of the .fam files generated in the previous step.
5. Genome-wide association analysis was performed using GEMMA v0.98.3:
   - Generate relatedness matrix:
   `gemma -bfile {famfile_prefix} -gk 1 -outdir {outdir} -o {outfile_prefix}`
   - Run linear mixed model:
   `gemma -bfile {famfile_prefix} -k {path_to_relatedness_matrix} -lmm 4 -outdir {outdir} -o {outfile_prefix}`

### 3. Gene ontology enrichment test (03\_)

#### Data

* **03_GO_fAstCal1.2-gene_table_Biomart_221108.txt**

  Gene table for _Astatotilapia calliptera_ (Ensembl Release 108, fAstCal1.2 genes) extracted from BioMart. Includes:
  * Gene stable ID
  * Transcript stable ID
  * Gene name
  * Transcription start site (TSS)
  * Gene start (bp)
  * Gene end (bp)
  * NCBI gene ID
  * Chromosome/scaffold name
* **03_GO_genes_GWAS_callset_ENSEMBL-genenames.txt**

  List of all genes annotated in the VCF used in the GWAS for relative eye size; serves as the gene universe for enrichment testing.
* **03_GO_genes_top001pct_25kbrange_ENSEMBL-genenames.txt**

  List of genes with TSS within 25 kb of top 0.01% GWAS outlier SNPs.

#### Notebooks

* **03_GO_enrichment_analysis_topGO.ipynb**

  Gene ontology enrichment test using the R package topGO. Takes as input the gene lists provided.

#### Workflow

1. The GWAS callset was annotated using:
   `snpEff -Xmx16g fAstCal1.2.99 {VCF} > {VCF_ann}`
2. All genes in the annotated GWAS callset were extracted using:
   `SnpSift extractFields -s ',' {VCF_ann} "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].BIOTYPE" > outfile.txt`
   The unique gene IDs from the first column were saved as "03_GO_genes_GWAS_callset_ENSEMBL-genenames.txt" and used as _gene universe_ for the enrichment test.
4. Using the BioMart gene table provided, a list of _test genes_ was generated, including genes with TSS within 25 kb of GWA outlier SNPs. This list is available as "03_GO_genes_top001pct_25kbrange_ENSEMBL-genenames.txt". Genomic position of GWA outlier SNPs are provided in Supplementary Data 5. 

### 4. Overrepresentation of genomic annotations (04\_) and selection tests (05\_)

#### Scripts

* **04_genomic_annotations_MAF-matched_permutations.py**

  Performs an overrepresentation test for genomic annotations in GWAS outlier SNPs using MAF-matched permutations.
* **05_rehh_xpehh_slurm.R**

  Computes the cross-population EHH statistic (XP-EHH) between two populations.
* **05_rehh_ihs_slurm.R**

  Computes the integrated haplotype score statistic (iHS) for a population.
* **05_rehh_xpehh_jobsubmission.slm**

  Example of SLURM submission script for XP-EHH and iHS calculation using the scripts above. Processes chromosomes in parallel.
* **05_xpehh_get_random_genome-wide-distr_max.py**

  Generates an empirical distribution of extreme XP-EHH scores by sampling random windows across the genome. It takes the output of "05_rehh_xpehh_slurm.R" as input.

#### Workflow
In preparation for running the selection test scripts (05_\*), the VCF files (**00_VCF/**) need to be subset for samples of the two populations tested: _D_. 'bigeye black dorsal' and _D. limnothrissa_. Sample IDs of these species can be found in Supplementary Data 3 under 'sequence_id'.
  
  Note: The scripts provided are formatted to be run by SLURM and take arguments from the submission script (**\*.slm**). The R scripts can be run this way by a SLURM-managed cluster or alternatively: 1) R scripts can be modified to include definitions of file paths and parameters, 2) Variables can be specifed as command-line arguments (e.g., `Rscript 05_rehh_xpehh_slurm.R <pop1> <pop2> <outdir> <chrom> <vcf_file_pop1> <vcf_file_pop2>`).

### 5. Gene expression analysis (06\_)

#### Data

* **06_fAstCal1.2/**

  Directory containing the reference FASTA and genomic annotation (GTF) for *Astatotilapia calliptera* (NCBI RefSeq assembly GCF_900246225.1), used for RNAseq read alignment. Retrieved on 27-10-2023.
* **06_rnaseq_count-matrix_htseq-count.txt**

  Tab-delimited gene count matrix used as input for differential expression analysis.
* **06_rnaseq_DESeq2_metadata.txt**

  Sample metadata used in the differential expression analysis (06_rnaseq_6_DEA_DESeq2.R)
* **06_fAstCal1.2-opsins.txt**

  Table with _A. calliptera_ opsin gene information, including:
  * Gene symbol
  * Gene ID
  * Genomic location
  * Start (bp)
  * End (bp)
* Raw RNAseq reads are available at [DOI:].

#### Scripts
* **06_rnaseq_1_fastqc.slm**

  Runs FastQC for quality control of raw RNA-seq reads; outputs individual QC reports per sample.
* **06_rnaseq_2_STAR_genome_index.slm**


  Generates STAR genome index from the reference FASTA and GTF files.
* **06_rnaseq_3_STAR_alignment.slm**

  Aligns raw sequencing reads to the reference genome using STAR.
* **06_rnaseq_4_htseq-count.slm**

  Produces gene count files per sample from aligned reads using htseq-count.
* **06_rnaseq_5_count-matrix-from-htseq-count-output.sh**

  Merges htseq-count outputs into a single count matrix for downstream analysis in DESeq2.
* **06_rnaseq_6_DEA_DESeq2.R**

  Performs differential expression analysis using the R package DESeq2. Input: "06_rnaseq_count-matrix_htseq-count.txt"

  Note: files with a .slm extension are Slurm submission scripts used to run each step of the RNA-seq pipeline on a high-performance computing (HPC) cluster. Relevant commands can be extracted and adapted from them.
 
#### Notebooks
* **06_rnaseq_opsin_expression.ipynb**

  Calculates relative opsin gene expression (in transcript per million, TPM)

#### Workflow
The scripts are provided as templates that can be run sequentially to generate a gene count matrix. All scripts reference files and directories relative to a base directory defined by `path_to_dir` with the following structure:
path_to_dir/
├── metadata/             # See Supplementary Data 4
├── raw_data/             # Raw FASTQ files (available from XXX)
├── fastqc/               # Output directory for FastQC results
├── STAR/
│   ├── 06_fAstCal1.2/    # Reference genome and GTF
│   ├── genome_index/     # STAR index files
│   └── alignment/        # STAR alignment output
└── HTSeq/
    ├── counts/           # Final gene count outputs (per sample and merged matrix)
    
    
  The final output is a tab-delimited gene count matrix generated by htseq-count, summarizing the number of reads mapping to coding sequences (CDS) per gene. Rows correspond to gene IDs (from the 'gene_id' attribute in the GTF file) and columns to individual samples. This file is used as input for the differential gene expression analysis ("06_rnaseq_6_DEA_DESeq2.R") and the calculation of relative opsin expression (notebook "06_rnaseq_opsin_expression").
