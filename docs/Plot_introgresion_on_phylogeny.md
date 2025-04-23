Maize Phylogeny Analysis with Multiple Variants
================
Maize Genetics Lab
2025-04-23

- [1. Setting up the Environment](#1-setting-up-the-environment)
- [2. Reading Input Files](#2-reading-input-files)
- [3. Creating the Name Mapping
  Table](#3-creating-the-name-mapping-table)
- [4. Trimming the Alignment](#4-trimming-the-alignment)
- [5. Extracting Variant Information](#5-extracting-variant-information)
- [6. Building a Phylogenetic Tree](#6-building-a-phylogenetic-tree)
- [7. Visualizing the Full Tree with
  ggtree](#7-visualizing-the-full-tree-with-ggtree)
- [8. Creating a Donor Subset
  Analysis](#8-creating-a-donor-subset-analysis)
- [9. Project File Structure](#9-project-file-structure)
- [10. Advantages of ggtree Over Base R
  Plotting](#10-advantages-of-ggtree-over-base-r-plotting)
- [11. Understanding the Key
  Functions](#11-understanding-the-key-functions)
- [12. Conclusion](#12-conclusion)
- [13. Troubleshooting Tips](#13-troubleshooting-tips)

This tutorial walks through a complete workflow for processing maize
genetic data, focusing on:

1.  Relabeling sequence identifiers in a FASTA alignment with meaningful
    taxonomic labels
2.  Extracting genotype information at multiple variant positions (I211V
    and A204T)
3.  Creating a phylogenetic tree with ancestry and variant information
4.  Rotating the tree to position the B73 reference genome at the top

## 1. Setting up the Environment

First, let’s install and load the necessary packages:

``` r
# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

# Load packages
library(Biostrings)  # For FASTA manipulation
library(tidyverse)   # For data manipulation

# Load custom functions from targseq package
# This package contains our specialized functions for maize analysis
library(targseq)
```

## 2. Reading Input Files

Before we begin, let’s understand the input files we’ll be working with:

    ## extdata/

    ## ├── ancestry_mismatch.csv            # Input: Ancestry call mismatches with tree position per sample (if any)

    ## ├── donor_data.csv                   # Input: Environmental data from donors

    ## ├── hpc1_aligned.fasta               # Input: Multiple sequence alignment of the hpc1 gene

    ## │                                    # Contains DNA sequences from multiple maize accessions

    ## │                                    # Technical IDs as sequence headers

    ## │                                    # ~3000 bp per sequence

    ## │                                    # Aligned using MAFFT

    ## │

    ## └── seqid_label.csv                 # Input: Metadata mapping file

    ##                                      # Links technical sequence IDs to meaningful labels

    ##                                      # Contains columns: seqid, label_1, label_2, label_3

    ##                                      # Includes ancestry information (ancestry_call)

    ##                                      # Classifies samples as 'Donor' or 'Recurrent'

Now, let’s set up our directory structure and read our input files:

``` r
# Create a project directory if it doesn't exist
project_dir <- "tree_plotting"
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}

# Set file paths within the project directory
aln_file <-  system.file("extdata", "hpc1_aligned.fasta", package="targseq")

metadata_file <- system.file("extdata", "seqid_label.csv", package="targseq")


# Read the mapping table with taxonomic metadata
metadata <- read.csv(metadata_file)

# Display the first few rows to understand the structure
head(metadata)

# Read the FASTA alignment using Biostrings
original_alignment <- readDNAStringSet(aln_file)

# View information about the alignment
cat("Number of sequences:", length(original_alignment), "\n")
cat("Sequence length:", width(original_alignment)[1], "bp\n")

# Display the first few sequence names
head(names(original_alignment))
```

## 3. Creating the Name Mapping Table

Now we’ll prepare our sequence names and map them to meaningful labels:

``` r
# Extract sequence IDs from the alignment and remove description fields
seq_ids <- names(original_alignment)
trimmed_ids <- gsub("\\s.*", "", seq_ids, perl=TRUE)

# Verify sequence ID counts match
cat("Original sequence count:", length(original_alignment), "\n")
cat("Trimmed IDs count:", length(trimmed_ids), "\n")

# Assign the trimmed IDs to the alignment object
names(original_alignment) <- trimmed_ids

# Create a mapping dataframe that connects sequence IDs with taxonomic labels
# We also keep track of the original order in the alignment
name_map <- data.frame(seqid = trimmed_ids) %>%
  mutate(aln_order = row_number()) %>%
  inner_join(metadata)

# Filter to include only sequences that exist in our alignment
name_map <- name_map %>%
  dplyr::filter(seqid %in% names(original_alignment))
rownames(name_map) <- name_map$label_3

# Report on mapping success
cat("Sequences successfully mapped:", nrow(name_map), "\n")

# Filter the alignment to keep only sequences we have metadata for
filtered_alignment <- original_alignment[name_map$seqid]

# Now replace the technical IDs with meaningful labels from label_3 column
names(filtered_alignment) <- name_map$label_3

# Store counts for later use
n_taxa <- length(filtered_alignment)
n_pos <- width(filtered_alignment)[1]

# Examine our newly labeled sequences
head(names(filtered_alignment))
```

## 4. Trimming the Alignment

Next, we’ll trim the alignment to focus on the region of interest:

``` r
# Load GenomicRanges for efficient sequence manipulation
library(GenomicRanges)

# Define trimming ranges (starting from position 2226 to the end)
trimmedRanges <- GRanges(
  seqnames = names(filtered_alignment),
  ranges = IRanges(
    start = rep(2226, n_taxa),
    end = rep(n_pos, n_taxa)
  )
)

# Perform the trimming operation
trimmed_alignment <- filtered_alignment[trimmedRanges]

# Verify the trimmed alignment
cat("Trimmed alignment length:", width(trimmed_alignment)[1], "bp\n")
head(names(trimmed_alignment))

# Writing the Processed Alignment to a New File
# Define output file path
output_alignment <- file.path(project_dir, "hpc1_nice_labels.fasta")

# Write the renamed alignment to a new file
writeXStringSet(trimmed_alignment, output_alignment, format="fasta")
cat("Renamed alignment written to:", output_alignment, "\n")
```

## 5. Extracting Variant Information

Now we can extract information for both variants using our specialized
function from the targseq package:

``` r
# Define our variants of interest with their patterns
variants <- data.frame(
  name = c("A204T","I211V"),
  pattern = c("GCCGTGGCGTGGCGC", "ATCACCCGC")
)

# Extract variant genotypes using our custom function
# get_variant_gt is a function from the targseq package that identifies 
# genotype at specific positions based on surrounding sequence patterns
variant_data <- get_variant_gt(variants, trimmed_alignment)

# Convert variant calls to standardized REF/ALT format
# This makes visualization more consistent
variant_data$A204T <- c("ALT", "REF")[as.factor(variant_data$A204T)]
variant_data$I211V <- c("REF", "ALT")[as.factor(variant_data$I211V)]

# Ensure row names match the sequence names for proper mapping
rownames(variant_data) <- names(trimmed_alignment)
```

## 6. Building a Phylogenetic Tree

Now we’ll build and visualize a phylogenetic tree from our processed
alignment:

``` r
# Load packages for phylogenetic analysis
library(ape)        # Basic phylogenetics functions
library(phangorn)   # Advanced phylogenetics 
library(phytools)   # Additional tree utilities

# Read the renamed alignment for phylogenetic analysis
hpc1_aln <- read.dna(output_alignment, format="fasta")

# Check labels in the alignment
labels(hpc1_aln)

# Convert to phangorn's phyDat format for phylogenetic analysis
hpc1_phyDat <- phyDat(hpc1_aln, type = "DNA", levels = NULL)

# Calculate distance matrix using JC69 model
# Note: You could use modelTest() to find the best model for your data
# mt <- modelTest(hpc1_phyDat)
dna_dist <- dist.ml(hpc1_phyDat, model="JC69")

# Build UPGMA tree and ladderize it for better visualization
hpc1_UPGMA <- ladderize(upgma(dna_dist), right = FALSE)

# Simple plot to check the tree
plot(hpc1_UPGMA, cex = 0.7, main = "UPGMA Tree Before Rotation")
```

## 7. Visualizing the Full Tree with ggtree

Now we’ll create a comprehensive visualization of our phylogenetic tree
with variant and ancestry information using ggtree, plus add
environmental data as a heatmap.

``` r
# Load required packages for visualization
library(ggtree)     # For tree visualization
library(ggplot2)    # For plotting
library(ggnewscale) # For multiple color scales in the same plot
library(colorspace)
library(GenomicRanges) # For handling genomic coordinates
```

Let’s build our visualization step by step:

``` r
# Rotate the tree to put B73 at the top using the pivot_on function
rotated_tree <- pivot_on(hpc1_UPGMA, "Zm_B73")

# Prepare ancestry information for the tree
ancestry <- name_map[, c("label_3", "ancestry_call"), drop=FALSE]
rownames(ancestry) <- name_map$label_3

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

# Display the basic tree with ancestry information
plot(p1)
```

Now let’s add environmental data as a heatmap alongside the tree:

``` r
# Create a dataframe with sequence IDs in their original order
alignment_order <- data.frame(seqid=names(trimmed_alignment))

name_map$seqid <- name_map$label_3
# Load  environmental data

# Prepare environmental parameters data
# Here we're using elevation data from our metadata

donor_file <-  system.file("extdata", "donor_data.csv", package="targseq")

donor_data <- read.csv(donor_file ) 

# Join all metadata while preserving the original alignment order
taxa_info  <- alignment_order %>%
  left_join(name_map) %>%           
  left_join(donor_data, by=c(donor_accession= "donor_id")) %>%
  select(seqid,donor_accession, label=label_3, ancestry_call, elevation)

env_par <- taxa_info [, c("elevation"), drop=FALSE]
colnames(env_par) <- c("Elevation")
rownames(env_par) <- taxa_info$seqid

# Add elevation data as a heatmap next to the tree
p2 <- gheatmap(p1, 
               data = env_par[, "Elevation", drop=FALSE], 
               offset = 0.025, width = 0.15, 
               colnames = TRUE, colnames_angle = 45, 
               colnames_offset_y = 0.4, hjust = 0,
               colnames_position = "top") +
  # Use diverging color palette centered at 1500m
  scale_fill_continuous_divergingx(name = "Elevation (masl)", 
                                   palette = 'RdBu', 
                                   mid = 1500, 
                                   n_interp = 25) +
  theme(legend.position = c(0.25, 0.5))

# Display the tree with elevation heatmap
plot(p2)
```

Next, let’s add our variant data as another heatmap layer:

``` r
# Define color palette for REF/ALT states
ref_alt_colors <- c("REF" = "tomato", "ALT" = "royalblue")

# Add variant data as another heatmap
# We use new_scale_fill() to create a separate color scale from the elevation data
p3 <- gheatmap(p2 + new_scale_fill(), 
               variant_data, 
               offset = 0.045, width = 0.12, 
               colnames = TRUE, colnames_angle = 45, 
               colnames_offset_y = 0.65, hjust = 0,
               colnames_position = "top") +
  scale_fill_manual(values = ref_alt_colors, name = "Allele") +
  theme(legend.position = c(0.25, 0.5))

# Display the tree with both elevation and variant heatmaps
plot(p3)
```

Finally, let’s add a sequence alignment view to complete our
visualization:

``` r
# Trim the alignment to focus on the region of interest
n_taxa <- length(trimmed_alignment) 
n_pos <- width(trimmed_alignment)[1]

# Make sure our alignment is ordered to match the tree
ordered_alignment <- trimmed_alignment[rotated_tree$tip.label]

# Convert to DNAbin format for msaplot
bin_alignment <- as.DNAbin(ordered_alignment)

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

# Display the complete visualization
plot(alignment_tree)

# Save the visualization as a PNG file
png_file <- file.path(project_dir, "target_sequencing_tree_alignment.png")
ggsave(alignment_tree, file = png_file, height = 12, width = 8)
cat("Full tree visualization complete. Plot saved as:", png_file, "\n")
```

## 8. Creating a Donor Subset Analysis

Often we want to focus our analysis on a specific subset of taxa. In
this case, we’ll create a subset that includes only taxa with donor
ancestry plus the B73 reference genome.

``` r
# First, identify B73 in our dataset
# Look for "B73" in the label column
b73_label <- grep("B73", name_map$label_3, value = TRUE)[1]
cat("B73 reference identified as:", b73_label, "\n")

# Ancestry miscalled lines
seqid_label_file <-  system.file("extdata", "seqid_label.csv", package="targseq")

seqid_label<- read.csv(seqid_label_file) 

ancestry_mismatch_file <-  system.file("extdata", "ancestry_mismatch.csv", package="targseq")

                      
ancestry_miscall <- read.csv(ancestry_mismatch_file) %>% 
  filter(locus=="hpc1") %>% 
  dplyr::select(fastq_sample ) %>%
  inner_join(seqid_label)

# Create a subset of taxa with donor ancestry plus B73
donor_subset <- name_map %>%
  mutate(label=label_3) %>%
  filter(ancestry_call == "Donor" | label == b73_label) %>%
  # filter ancestry_miscalls
  filter(!seqid %in% ancestry_miscall$seqid)

# Check how many taxa we have in our subset
cat("Number of taxa in donor subset:", nrow(donor_subset), "\n")
```

Now that we have our subset defined, let’s extract the corresponding
sequences and build a new tree:

``` r
# Extract the subset of sequences from our alignment
donor_alignment <- trimmed_alignment[donor_subset$label_3]

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
```

Let’s visualize the donor subset tree with the same approach we used for
the full tree:

``` r
# Extract variant data for the donor subset
donor_variant_data <- variant_data[donor_subset$label_3, ]

# Extract environmental data for the donor subset
donor_env_data <- env_par[donor_subset$label_3, , drop=FALSE]

# Create ancestry information for the donor subset
donor_ancestry <- data.frame(
  label_3 = donor_subset$label_3,
  ancestry_call = donor_subset$ancestry_call
)
rownames(donor_ancestry) <- donor_subset$label_3

# Create the base tree with ancestry information
donor_p1 <- ggtree(donor_tree, ladderize = FALSE) %<+% donor_ancestry +
  # Add colored points indicating ancestry
  geom_tippoint(aes(color = ancestry_call), position = position_nudge(x = 0.001)) +
  # Add tip labels
  geom_tiplab(size = 3, offset = 0.002, align = TRUE, linesize = 0.1) +
  # Add vertical expansion for better spacing
  ggtree::vexpand(.1, 1) +
  # Use the custom ancestry color palette
  scale_color_manual(name="Ancestry", values = ancestry_palette) +
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
  scale_fill_continuous_divergingx(name = "Elevation (masl)", 
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
  theme(legend.position = c(0.25, 0.5))

# Create subset of donor alignment for msaplot
donor_bin_alignment <- as.DNAbin(donor_alignment)

# Add alignment view alongside the tree and heatmaps
donor_alignment_tree <- msaplot(donor_p3 + new_scale_fill(), 
                                fasta = donor_bin_alignment, 
                                offset = 0.06, 
                                width = 0.5) +
  guides(fill = "none") +
  theme(legend.position = c(0.25, 0.5))

# Display the donor subset visualization
plot(donor_alignment_tree)

# Save the donor subset visualization
donor_plot_file <- file.path(project_dir, "donor_subset_tree_alignment.png")
ggsave(donor_alignment_tree, file = donor_plot_file, height = 10, width = 9)
cat("Donor subset analysis complete. Plot saved as:", donor_plot_file, "\n")
```

The donor subset analysis allows us to focus specifically on the
relationships between B73 and the donor accessions, making it easier to
identify patterns of introgression or shared ancestry. By highlighting
B73 in gold, we can quickly identify it as our reference genome in the
visualization.

This subset approach is particularly useful when working with large
datasets where the full tree might be too crowded to interpret easily.
It also allows us to examine how specific groups of taxa (in this case,
donor accessions) relate to a reference genome of interest.

## 9. Project File Structure

Let’s examine the complete file structure with all our inputs and
outputs:

    ## tree_plotting/

    ## ├── hpc1_aligned.fasta              # Input: Original alignment file

    ## ├── seqid_label.csv                 # Input: Metadata mapping file

    ## ├── hpc1_nice_labels.fasta          # Output: Alignment with renamed sequences

    ## ├── hpc1_donor_subset.fasta         # Output: Subset alignment with donor taxa + B73

    ## ├── hpc1_UPGMA.tre                  # Output: Full phylogenetic tree in Newick format

    ## ├── hpc1_donor_subset.tre           # Output: Donor subset tree in Newick format

    ## ├── target_sequencing_tree_alignment.png # Output: Full tree with visualization

    ## └── donor_subset_tree_alignment.png      # Output: Donor subset visualization

## 10. Advantages of ggtree Over Base R Plotting

The ggtree package offers several advantages over base R plotting
functions for phylogenetic trees:

1.  **Grammar of Graphics Approach**: Uses the familiar ggplot2 syntax
    for consistent and extensible visualizations.

2.  **Layered Visualizations**: Easily add multiple data layers
    (geom_tippoint, geom_hilight, etc.) to the tree.

3.  **Multiple Tree Layouts**: Supports rectangular, circular, fan, and
    other tree layouts without rewriting code.

4.  **Integration with Other Data**: The `%<+%` operator makes it easy
    to map external data onto the tree.

5.  **Publication-Ready Figures**: Creates high-quality vector graphics
    with customizable aesthetics.

6.  **Highlighting and Annotation**: Easily highlight specific clades or
    nodes of interest.

7.  **Multiple Output Formats**: Compatible with all ggplot2 output
    formats (PDF, PNG, SVG, etc.).

8.  **Heatmap Integration**: The gheatmap function allows easy addition
    of heatmaps alongside trees.

## 11. Understanding the Key Functions

Let’s briefly review the key functions from our package that were used
in this workflow:

### 1. Variant Detection Function

The `get_variant_gt()` function identifies variant positions in DNA
sequences by: - Scanning for specific sequence patterns around known
variant sites - Extracting the nucleotide at the variant position -
Converting the calls to a standardized format (REF/ALT)

### 2. Tree Rotation Function

The `pivot_on()` function rotates a phylogenetic tree to position a
specific taxon at the top by: - Finding the node corresponding to the
specified taxon - Traversing the tree structure to determine necessary
rotations - Preserving the original tree topology and branch lengths

### 3. Tree Visualization Functions

The ggtree package provides multiple functions for tree visualization: -
`ggtree()`: Creates the basic tree structure - `geom_tippoint()`: Adds
points at the tips representing data values - `geom_tiplab()`: Adds
labels at the tips - `gheatmap()`: Adds heatmaps alongside the tree -
`msaplot()`: Adds sequence alignment views

## 12. Conclusion

This tutorial demonstrated how to use our specialized functions with the
ggtree package to create advanced visualizations of phylogenetic trees
with variant and ancestry information. We created:

1.  Full tree visualizations with nucleotide and REF/ALT variant coding
2.  Subset analyses focusing on donor taxa plus the B73 reference
3.  Integrated heatmaps for variant and environmental data visualization
4.  Sequence alignment views that provide direct visualization of the
    underlying data

These visualizations help reveal patterns of ancestry and genetic
variation that would be difficult to detect through other means.

## 13. Troubleshooting Tips

- **Missing data issues**: Use the `%<+%` operator from ggtree to
  properly map data to the tree tips.
- **Layout problems**: If the tree is too crowded, try adjusting the
  `width` and `height` parameters in `ggsave()`.
- **Label overlaps**: Use `geom_tiplab2()` with an appropriate `offset`
  value to prevent overlapping labels.
- **Color issues**: Define custom color palettes using
  `scale_color_manual()` or `scale_fill_manual()`.
- **Legend problems**: Adjust the legend position and box orientation
  with `theme()` options.
