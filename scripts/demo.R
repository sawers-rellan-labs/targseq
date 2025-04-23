#' Maize Phylogeny Analysis with Variant and Environmental Data
#' 
#' This script demonstrates how to analyze maize genetic data by:
#' 1. Reading sequence alignment and metadata
#' 2. Extracting variant information
#' 3. Building a phylogenetic tree
#' 4. Visualizing the tree with variant states, ancestry, and environmental data
#' 5. Creating a subset analysis focusing on donor ancestry taxa
#' 
#' Authors: Maize Genetics Lab
#' Version: 1.1

# Load required libraries
library(targseq)      # Custom package for targeted sequencing analysis
library(Biostrings)   # For handling biological sequences
library(phangorn)     # For phylogenetic analysis
library(ape)          # For phylogenetic tree operations
library(ggtree)       # For tree visualization
library(ggplot2)      # For plotting
library(colorspace)   # For color manipulation
library(dplyr)        # For data manipulation
library(ggnewscale)   # For multiple scales in ggplot
library(GenomicRanges) # For handling genomic coordinates

#===============================================================================
# 1. Read sequence alignment and metadata
#===============================================================================

# Create a directory for outputs if it doesn't exist
project_dir <- "tree_plotting"
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}

# Load sequence alignment from the targseq package
aln_file <- system.file("extdata", "hpc1_aligned.fasta", package="targseq")
alignment <- readDNAStringSet(aln_file)

# Create a dataframe with sequence IDs in their original order
alignment_order <- data.frame(seqid=names(alignment)) 

# Load metadata files from the targseq package
metadata_file <- system.file("extdata", "seqid_label.csv", package="targseq")
donor_file <- system.file("extdata", "donor_data.csv", package="targseq")
donor_data <- read.csv(donor_file)

# Read and join the metadata information
cat("Reading sample metadata...\n")
taxa_info <- read.csv(metadata_file)

# Join all metadata while preserving the original alignment order
taxa_info <- alignment_order %>%
  left_join(taxa_info) %>%           
  left_join(donor_data, by=c(donor_accession= "donor_id")) %>%
  select(seqid, label=label_3, ancestry_call, elevation)

# Rename sequences with meaningful labels
names(alignment) <- taxa_info$label

#===============================================================================
# 2. Extract variant genotypes
#===============================================================================

# Define variants to extract: A204T and I211V
# Each variant is identified by a unique DNA pattern
variants <- data.frame(
  name = c("A204T", "I211V"),
  pattern = c("GCCGTGGCGTGGCGC", "ATCACCCGC")
)

# Extract variant genotypes from the alignment
# (get_variant_gt is a function from the targseq package)
variant_data <- get_variant_gt(variants, alignment)

#===============================================================================
# 3. Recode variants as REF/ALT
#===============================================================================

# Recode A204T: A=REF, G=ALT
variant_data$A204T <- c("ALT", "REF")[as.factor(variant_data$A204T)]

# Recode I211V: G=REF, A=ALT
variant_data$I211V <- c("REF", "ALT")[as.factor(variant_data$I211V)]

# Set row names for easier integration with the tree
rownames(variant_data) <- taxa_info$label

#===============================================================================
# 4. Build phylogenetic tree
#===============================================================================

# Convert FASTA to phyDat format for phylogenetic analysis
phyDat <- phyDat(read.dna(aln_file, format="fasta"), type="DNA")
names(phyDat) <- taxa_info$label

# Calculate genetic distances using the JC69 model
dna_dist <- dist.ml(phyDat, model="JC69")

# Build a UPGMA tree and ladderize it
tree <- ladderize(upgma(dna_dist), right=FALSE)

# Rotate the tree to put B73 at the top
# pivot_on is a custom function that rotates the tree to place a specific taxon at the top
rotated_tree <- pivot_on(tree, "Zm_B73")

#===============================================================================
# 5. Create tree visualization with ggtree
#===============================================================================

# Add ancestry information to the tree
ancestry <- taxa_info[, c("label", "ancestry_call"), drop=FALSE]
rownames(ancestry) <- taxa_info$label

# Define color palette for ancestry
ancestry_palette <- c("Recurrent" = "tomato", "Donor" = "royalblue")

# Create the base tree with ancestry information
p1 <- ggtree(rotated_tree, ladderize = FALSE) %<+% ancestry +
  # Add colored points indicating ancestry
  geom_tippoint(aes(color = ancestry_call), position = position_nudge(x = 0.0015)) +
  # Add tip labels
  geom_tiplab(size = 2.5, offset = 0.002, align = TRUE, linesize = 0.1) +
  # Add vertical expansion for better spacing
  ggtree::vexpand(.1, 1) +
  # Use the custom ancestry color palette
  scale_color_manual(name="Ancestry call", values = ancestry_palette) +
  # Position the legend
  theme(legend.position = c(0.25, 0.5)) +
  # Use tree theme for clean visualization
  theme_tree2()

#===============================================================================
# 6. Add environmental data as a heatmap
#===============================================================================

# Prepare environmental parameters data
env_par <- taxa_info[, c("elevation"), drop=FALSE]
colnames(env_par) <- c("Elevation")
rownames(env_par) <- taxa_info$label

# Add elevation data as a heatmap next to the tree
p2 <- gheatmap(p1, 
               data = env_par[, "Elevation", drop=FALSE], 
               offset = 0.025, width = 0.15, 
               colnames = TRUE, colnames_angle = 45, 
               colnames_offset_y = 0.4, hjust = 0,
               colnames_position = "top") +
  # Use diverging color palette centered at 1500m
  scale_fill_continuous_divergingx(name = "Donor Elevation (masl)", 
                                   palette = 'RdBu', 
                                   mid = 1500, 
                                   n_interp = 25) +
  theme(legend.position = c(0.25, 0.5))

#===============================================================================
# 7. Add variant data as another heatmap
#===============================================================================

# Define color palette for REF/ALT states
ref_alt_colors <- c("REF" = "tomato", "ALT" = "royalblue")

# Add variant data as another heatmap
# Use new_scale_fill() to create a separate color scale from the elevation data
p3 <- gheatmap(p2 + new_scale_fill(), 
               variant_data, 
               offset = 0.045, width = 0.12, 
               colnames = TRUE, colnames_angle = 45, 
               colnames_offset_y = 0.65, hjust = 0,
               colnames_position = "top") +
  scale_fill_manual(values = ref_alt_colors, name = "Allele") +
  theme(legend.position = c(0.25, 0.5))

#===============================================================================
# 8. Add sequence alignment view
#===============================================================================

# Trim the alignment to focus on the region of interest
n_taxa <- length(alignment) 
n_pos <- width(alignment)[1]

# Create genomic ranges for trimming (starting from position 2226)
trimmedRanges <- GRanges(
  seqnames = names(alignment),
  ranges = IRanges(
    start = rep(2226, n_taxa),
    end = rep(n_pos, n_taxa)
  )
)

# Perform the trimming operation
trimmed_alignment <- alignment[trimmedRanges]

# Order the alignment to match the tree
trimmed_alignment <- trimmed_alignment[rotated_tree$tip.label]

# Convert to DNAbin format for msaplot
bin_alignment <- as.DNAbin(trimmed_alignment)

# Verify alignment order matches the tree
all(names(bin_alignment) == rotated_tree$tip.label)

# Add alignment view alongside the tree and heatmaps
# Use new_scale_fill() again to create a separate color scale
alignment_tree <- msaplot(p3 + new_scale_fill(), 
                          fasta = bin_alignment, 
                          offset = 0.06, 
                          width = 0.5) +
  guides(fill = "none") +
  theme(legend.position = c(0.25, 0.5))

#===============================================================================
# 9. Save the final plot
#===============================================================================

# Save the visualization as a PNG file
png_file <- file.path(project_dir,"target_sequencing_tree_alignment.png")
ggsave(alignment_tree, file = png_file, height = 12, width = 8)

# Print a message confirming the analysis is complete
cat("Full analysis complete. Plot saved as 'maize_tree_alignment.png'\n")

#===============================================================================
# 10. Create and analyze subset of taxa with donor ancestry plus B73
#===============================================================================



# First, identify B73 in our dataset
# Look for "B73" in the label column
b73_label <- grep("B73", taxa_info$label, value = TRUE)[1]
cat("B73 reference identified as:", b73_label, "\n")

# Ancestry miscalled lines
seqid_label_file <-  system.file("extdata", "seqid_label.csv", package="targseq")

seqid_label<- read.csv(seqid_label_file) 

ancestry_mismatch_file <-  system.file("extdata", "ancestry_mismatch.csv", package="targseq")

ancestry_miscall <- read.csv(ancestry_mismatch_file) %>% 
  filter(locus=="hpc1") %>%
  inner_join(seqid_label, by =c(locus="gene","fastq_sample"))

ancestry_miscall

# Create a subset of taxa with donor ancestry plus B73
donor_subset <- taxa_info %>%
  filter(ancestry_call == "Donor" | label == b73_label) %>%
  # filter ancestry_miscalls
  filter(!seqid %in% ancestry_miscall$seqid)


# Check how many taxa we have in our subset
cat("Number of taxa in donor subset:", nrow(donor_subset), "\n")

# Extract the subset of sequences from our alignment
donor_alignment <- trimmed_alignment[donor_subset$label]

# Write the subset alignment to a file
donor_alignment_file <- file.path(project_dir, "hpc1_donor_subset.fasta")
writeXStringSet(donor_alignment, donor_alignment_file, format="fasta")
cat("Donor subset alignment written to:", donor_alignment_file, "\n")

# Build a new tree for the donor subset
donor_aln <- read.dna(donor_alignment_file, format="fasta")
donor_phyDat <- phyDat(donor_aln, type = "DNA", levels = NULL)
donor_dist <- dist.ml(donor_phyDat, model="JC69")
donor_UPGMA <- ladderize(upgma(donor_dist), right = FALSE)

# Rotate the donor tree to put B73 at the top
donor_tree <- pivot_on(donor_UPGMA, b73_label)

# Save the donor subset tree
donor_tree_file <- file.path(project_dir, "hpc1_donor_subset.tre")
write.tree(donor_tree, file = donor_tree_file)
cat("Donor subset tree saved to:", donor_tree_file, "\n")

#===============================================================================
# 11. Create visualizations for the donor subset
#===============================================================================

# Extract variant data for the donor subset
donor_variant_data <- variant_data[donor_subset$label, ]

# Extract environmental data for the donor subset
donor_env_data <- env_par[donor_subset$label, , drop=FALSE]

# Create ancestry information for the donor subset
donor_ancestry <- data.frame(
  label = donor_subset$label,
  ancestry_call = donor_subset$ancestry_call
)
rownames(donor_ancestry) <- donor_subset$label

# Create the base tree with ancestry information
donor_p1 <- ggtree(donor_tree, ladderize = FALSE) %<+% donor_ancestry +
  # Add colored points indicating ancestry
  geom_tippoint(aes(color = ancestry_call), position = position_nudge(x = 0.001)) +
  # Add tip labels
  geom_tiplab(size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
  # Add vertical expansion for better spacing
  ggtree::vexpand(.1, 1) +
  # Use the custom ancestry color palette
  scale_color_manual(name="Ancestry call", values = ancestry_palette) +
  # Highlight B73 reference
  geom_hilight(node = which(donor_tree$tip.label == b73_label), 
               fill = "gold", alpha = 0.3) +
  # Position the legend
  theme(legend.position = c(0.25, 0.5)) +
  # Use tree theme for clean visualization
  theme_tree2() +
  # Add title
  ggtitle("Donor Ancestry Taxa with B73 Reference")

# Add environmental data as a heatmap
donor_p2 <- gheatmap(donor_p1, 
                     data = donor_env_data[, "Elevation", drop=FALSE], 
                     offset = 0.025, width = 0.15, 
                     colnames = TRUE, colnames_angle = 45, 
                     colnames_offset_y = 0.4, hjust = 0,
                     colnames_position = "top") +
  # Use diverging color palette centered at 1500m
  scale_fill_continuous_divergingx(name = "Donor Elevation (masl)", 
                                   palette = 'RdBu', 
                                   mid = 1500, 
                                   n_interp = 25) +
  theme(legend.position = c(0.25, 0.5))

# Add variant data as another heatmap
donor_p3 <- gheatmap(donor_p2 + new_scale_fill(), 
                     donor_variant_data, 
                     offset = 0.045, width = 0.12, 
                     colnames = TRUE, colnames_angle = 45, 
                     colnames_offset_y = 0.65, hjust = 0,
                     colnames_position = "top") +
  scale_fill_manual(values = ref_alt_colors, name = "Allele") +
  theme(legend.position = c(0.75, 0.5))

# Create subset of donor alignment for msaplot
donor_bin_alignment <- as.DNAbin(donor_alignment)

# Add alignment view alongside the tree and heatmaps
donor_alignment_tree <- msaplot(donor_p3 + new_scale_fill(), 
                                fasta = donor_bin_alignment, 
                                offset = 0.06, 
                                width = 0.5) +
  guides(fill = "none") +
  theme(legend.position = c(0.25, 0.75))

# Save the donor subset visualizationx
donor_plot_file <- file.path(project_dir, "donor_subset_tree_alignment.png")
ggsave(donor_alignment_tree, file = donor_plot_file, height = 10, width = 9)
cat("Donor subset analysis complete. Plot saved as:", donor_plot_file, "\n")



# Print a message confirming all analyses are complete
cat("\nAll analyses complete. Output files are in the  directory.\n")

