# Last modified: 14-01-2025

# This script calculates phylogenetic signal for relative eye size in Diplotaxodon species 
# using Pagel's lambda method from the phytools package.

# Input files required:
# - (1) tab-delimited file with species names and relative eye size measurements with columns:
#       - full_name: species (e.g., Diplotaxodon_limnothrissa)      
#       - ED.SL: relative eye size (eye diameter divided by standard length)
#
# This file can be derived from Supplementary Data 1 columns 'ED' (eye diameter), 'SL' (standard length), and 'full_name'.
#      
# - (2) Diplotaxodon species tree. Provided as "02_NJ_tree_diplotaxodon.txt"
#
# Below, (1) and (2) are read as 'eye_size_df' and 'diplo_tree', respectively.

#------------------------------------------------------------------------------
#                       1. Load required packages
#------------------------------------------------------------------------------
# Install phytools if not already installed
if (!requireNamespace("phytools", quietly = TRUE)) {
  install.packages("phytools")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Load packages
library(phytools)
library(dplyr)

# Check package version
packageVersion("phytools")

#------------------------------------------------------------------------------
#                       2. Import and prepare data
#------------------------------------------------------------------------------
# Load eye size data
eye_size_df <- read.csv("path_to_file/eyesize.txt",      # specify path accordingly
                        sep = '\t')

# Calculate mean eye size by species
mean_eye_size_sp <- aggregate(eye_size_df$ED.SL, 
                               by = list(species = eye_size_df$full_name), 
                               FUN = mean)
colnames(mean_eye_size_sp) <- c("species", "mean_eye_size")
mean_eye_size_sp$species <- as.character(mean_eye_size_sp$species)

# Load phylogenetic tree
diplo_tree <- read.newick(file = "path_to_file/02_NJ_tree_diplotaxodon.txt")    # specify path  

#------------------------------------------------------------------------------
#                    3. Match phenotype data with tree labels
#------------------------------------------------------------------------------
# Add '.0' suffix to species names to match tree labels
mean_eye_size_sp <- mean_eye_size_sp %>% 
  mutate(species = paste0(species, ".0"))

# The two Diplotaxodon 'holochromis' are not monophyletic. Define as separate species:
# Manually add Diplotaxodon holochromis.1 (cichlid6978788)
new_row <- data.frame(species = "Diplotaxodon_holochromis.1", mean_eye_size = 0.078769)
mean_eye_size_sp <- bind_rows(mean_eye_size_sp, new_row)

# Manually assign phenotype value to Diplotaxodon holochromis.0 (cichlid6978780)
mean_eye_size_sp$mean_eye_size[which(mean_eye_size_sp$species == "Diplotaxodon_holochromis.0")] <- 0.083108

# Get list of species with phenotype data
species_with_data <- unique(mean_eye_size_sp$species)

#------------------------------------------------------------------------------
#                          4. Prune the phylogenetic tree
#------------------------------------------------------------------------------
# Identify species in tree with no phenotype data available
list_sp_remove <- setdiff(diplo_tree$tip.label, species_with_data)

# Drop species with no eye size data from tree
diplo_tree_subset <- drop.tip(diplo_tree, list_sp_remove)

# D. similis not in the tree - remove
mean_eye_size_sp <- mean_eye_size_sp[which(!mean_eye_size_sp$species == 'Diplotaxodon_similis.0'), ]

#------------------------------------------------------------------------------
#                           5. Estimate phylogenetic signal
#------------------------------------------------------------------------------
# Create named vector of trait values for phylosig function
mean_eye_size_test <- setNames(mean_eye_size_sp$mean_eye_size, mean_eye_size_sp$species)

# Run phylogenetic signal test using Pagel's lambda
phylosig_lambda <- phylosig(tree = diplo_tree_subset, 
                           x = mean_eye_size_test,
                           method = "lambda")
# Print phylogenetic signal results
print(phylosig_lambda)

# Output results to a text file
#write.table(data.frame(lambda = phylosig_lambda$lambda,
#                       logL = phylosig_lambda$logL),
#            file = "phylogenetic_signal_results.txt", 
#            row.names = FALSE, quote = FALSE)