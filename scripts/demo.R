#!/usr/bin/env Rscript
#############################################################################
# targseq_demo.R - Demonstration of the targseq package
# 
# This script demonstrates how to use the targseq package for maize phylogenetic
# analysis with variant visualization.
#
# Author: Maize Genetics Lab
# Date: April 2025
#############################################################################

#===========================================================================
# 1. SETUP AND INSTALLATION
#===========================================================================

#---------------------------------------------------------------------------
# 1.1 Install required packages (if needed)
#---------------------------------------------------------------------------
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install targseq package from GitHub
# Only run this if you need to install or update the package
# Comment out after successful installation
# devtools::install_github("sawers-rellan-labs/targseq")

#---------------------------------------------------------------------------
# 1.2 Load required libraries
#---------------------------------------------------------------------------
# Main targseq package
library(targseq)

# Bioinformatics packages
library(Biostrings)  # For sequence handling
library(phangorn)    # For phylogenetic analysis
library(ape)         # For tree manipulation

# Visualization packages
library(ggtree)      # For tree visualization
library(ggplot2)     # For plotting
library(dplyr)       # For data manipulation

# Set output directory
project_dir <- "."
# Uncomment below to specify a different output directory
# project_dir <- "~/path/to/output"

# Check if directory exists and create if needed
if (!dir.exists(project_dir)) {
  dir.create(project_dir, recursive = TRUE)
}

#===========================================================================
# 2. DATA IMPORT & PREPROCESSING
#===========================================================================

#---------------------------------------------------------------------------
# 2.1 Load sequence alignment
#---------------------------------------------------------------------------
# Get example data path from package
aln_file <- system.file("extdata", "hpc1_aligned.fasta", package="targseq")

# Read the alignment
cat("Reading sequence alignment...\n")
alignment <- readDNAStringSet(aln_file)

# Display basic stats
cat(sprintf("Loaded %d sequences, each %d bp long\n", 
            length(alignment), width(alignment)[1]))

#---------------------------------------------------------------------------
# 2.2 Load metadata
#---------------------------------------------------------------------------
# Get metadata file path from package
metadata_file <- system.file("extdata", "seqid_label.csv", package="targseq")

# Read metadata 
cat("Reading sample metadata...\n")
taxa_info <- read.csv(metadata_file)

# Set rownames for easy indexing
rownames(taxa_info) <- taxa_info$seqid

# Rename ancestry column for consistency
taxa_info <- taxa_info %>%
  rename(ancestry_call = "founder_ancestry")

#===========================================================================
# 3. VARIANT DETECTION
#===========================================================================

#---------------------------------------------------------------------------
# 3.1 Define variants of interest
#---------------------------------------------------------------------------
cat("Detecting variants...\n")
# Define variants with their flanking sequence patterns
variants <- data.frame(
  name = c("A204T", "I211V"),
  pattern = c("GCCGTGGCGTGGCGC", "ATCACCCGC")
)

#---------------------------------------------------------------------------
# 3.2 Extract variant information
#---------------------------------------------------------------------------
# The get_variant_gt function finds these patterns in the sequences
# and extracts the variant nucleotide at the specified position
variant_data <- get_variant_gt(variants, alignment)

# Convert to standardized REF/ALT format for visualization
variant_data$A204T <- c("ALT", "REF")[as.factor(variant_data$A204T)]
variant_data$I211V <- c("REF", "ALT")[as.factor(variant_data$I211V)]

# Ensure row names match the sequence names
rownames(variant_data) <- names(alignment)

#===========================================================================
# 4. PHYLOGENETIC TREE CREATION
#===========================================================================

#---------------------------------------------------------------------------
# 4.1 Build phylogenetic tree
#---------------------------------------------------------------------------
cat("Building phylogenetic tree...\n")
# Convert to phyDat format for phylogenetic analysis
phyDat <- phyDat(data = read.dna(aln_file, format="fasta"), 
                 type = "DNA", levels = NULL)

# Calculate distance matrix using JC69 model
# Note: For more rigorous analysis, use modelTest() to find the best model
# mt <- modelTest(phyDat)
dna_dist <- dist.ml(phyDat, model="JC69")

# Build UPGMA tree and arrange branches for better visualization
tree <- ladderize(upgma(dna_dist), right = FALSE)

#---------------------------------------------------------------------------
# 4.2 Rotate tree to place B73 reference at top
#---------------------------------------------------------------------------
# The pivot_on function rotates the tree to place a specified taxon at the top
# In this case, we want B73 (reference genome) at the top
rotated_tree <- pivot_on(tree, "hpc1_B73")

#===========================================================================
# 5. VISUALIZATION
#===========================================================================

#---------------------------------------------------------------------------
# 5.1 Create basic tree visualization with variants
#---------------------------------------------------------------------------
cat("Creating visualizations...\n")
# Prepare taxa info for the specific tree
tree_taxa_info <- taxa_info[rotated_tree$tip.label, ]

# Create variant heatmap tree
tree_plot <- create_variant_heatmap_tree(
  tree = rotated_tree, 
  data = variant_data
)

# Define color palette for ancestry
ancestry_palette <- c("Recurrent" = "tomato", "Donor" = "royalblue")

#---------------------------------------------------------------------------
# 5.2 Add ancestry information to tree
#---------------------------------------------------------------------------
# Enhance tree with ancestry information using ggtree's %<+% operator
enhanced_tree_plot <- tree_plot %<+% tree_taxa_info +
  # Add colored points indicating ancestry
  geom_tippoint(aes(color = ancestry_call), 
                position = position_nudge(x = 0.0015)) +
  # Set color scheme for ancestry points
  scale_color_manual(values = ancestry_palette) +
  # Position the legend
  theme(legend.position = c(0.25, 0.5))

# Save the enhanced tree plot
output_file <- file.path(project_dir, "hpc1_tree_plot.pdf")
ggsave(enhanced_tree_plot, file = output_file,
       height = 12, width = 7, units = "in")
cat(sprintf("Saved tree plot to: %s\n", output_file))

#---------------------------------------------------------------------------
# 5.3 Create alternative visualization (tree with tip labels)
#---------------------------------------------------------------------------
# Create a simpler tree with just labels and ancestry points
simple_tree <- ggtree(rotated_tree, ladderize = FALSE) %<+% 
  tree_taxa_info + 
  geom_tiplab(size = 3) +
  geom_tippoint(aes(color = ancestry_call), 
                position = position_nudge(x = 0.0015)) +
  scale_color_manual(values = ancestry_palette)

#---------------------------------------------------------------------------
# 5.4 Add multiple sequence alignment to tree visualization
#---------------------------------------------------------------------------
# Create combined visualization with tree + alignment
alignment_tree <- msaplot(enhanced_tree_plot, 
                          fasta = aln_file, 
                          offset = 0.05,   # Space between tree and alignment
                          width = 0.6) +   # Width of alignment panel
  theme(legend.position = c(0.25, 0.5))

# Save the tree with alignment visualization
output_file <- file.path(project_dir, "hpc1_tree_alignment.png")
ggsave(alignment_tree, 
       file = output_file,
       height = 12, width = 7, units = "in")
cat(sprintf("Saved tree with alignment to: %s\n", output_file))

cat("Demo completed successfully!\n")