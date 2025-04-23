#' Maize Phylogeny Analysis
#' Version: 1.1

# Load libraries
library(targseq)
library(Biostrings)
library(phangorn)
library(ape)
library(ggtree)
library(ggplot2)
library(colorspace)
library(dplyr)
library(ggnewscale)
library(GenomicRanges)

# Create output directory
project_dir <- "tree_plotting"
if (!dir.exists(project_dir)) dir.create(project_dir)

# Read data
aln_file <- system.file("extdata", "hpc1_aligned.fasta", package="targseq")
alignment <- readDNAStringSet(aln_file)
alignment_order <- data.frame(seqid=names(alignment))

metadata_file <- system.file("extdata", "seqid_label.csv", package="targseq")
donor_file <- system.file("extdata", "donor_data.csv", package="targseq")
donor_data <- read.csv(donor_file)

# Process metadata
taxa_info <- read.csv(metadata_file) %>%
  right_join(alignment_order) %>%
  left_join(donor_data, by=c(donor_accession="donor_id")) %>%
  select(seqid, label=label_3, ancestry_call, elevation)

# Rename sequences
names(alignment) <- taxa_info$label

# Extract variant data
variants <- data.frame(
  name = c("A204T", "I211V"),
  pattern = c("GCCGTGGCGTGGCGC", "ATCACCCGC")
)
variant_data <- get_variant_gt(variants, alignment)

# Recode variants
variant_data$A204T <- c("ALT", "REF")[as.factor(variant_data$A204T)]
variant_data$I211V <- c("REF", "ALT")[as.factor(variant_data$I211V)]
rownames(variant_data) <- taxa_info$label

# Build tree
phyDat <- phyDat(read.dna(aln_file, format="fasta"), type="DNA")
names(phyDat) <- taxa_info$label
dna_dist <- dist.ml(phyDat, model="JC69")
tree <- ladderize(upgma(dna_dist), right=FALSE)
rotated_tree <- pivot_on(tree, "Zm_B73")

# Create tree visualization
ancestry_palette <- c("Recurrent" = "tomato", "Donor" = "royalblue")
ancestry <- taxa_info[, c("label", "ancestry_call"), drop=FALSE]
rownames(ancestry) <- taxa_info$label

p1 <- ggtree(rotated_tree, ladderize = FALSE) %<+% ancestry +
  geom_tippoint(aes(color = ancestry_call), position = position_nudge(x = 0.0015)) +
  geom_tiplab(size = 2.5, offset = 0.002, align = TRUE, linesize = 0.1) +
  ggtree::vexpand(.1, 1) +
  scale_color_manual(name="Ancestry call", values = ancestry_palette) +
  theme(legend.position = c(0.25, 0.5)) +
  theme_tree2()

# Add environmental data
env_par <- taxa_info[, c("elevation"), drop=FALSE]
colnames(env_par) <- c("Elevation")
rownames(env_par) <- taxa_info$label

p2 <- gheatmap(p1, env_par[, "Elevation", drop=FALSE], 
               offset = 0.025, width = 0.15, 
               colnames = TRUE, colnames_angle = 45,
               colnames_offset_y = 0.4, hjust = 0,
               colnames_position = "top") +
  scale_fill_continuous_divergingx(name = "Elevation (masl)", 
                                   palette = 'RdBu', 
                                   mid = 1500, 
                                   n_interp = 25) +
  theme(legend.position = c(0.25, 0.5))

# Add variant data
ref_alt_colors <- c("REF" = "tomato", "ALT" = "royalblue")
p3 <- gheatmap(p2 + new_scale_fill(), variant_data, 
               offset = 0.045, width = 0.12, 
               colnames = TRUE, colnames_angle = 45, 
               colnames_offset_y = 0.65, hjust = 0,
               colnames_position = "top") +
  scale_fill_manual(values = ref_alt_colors, name = "Allele") +
  theme(legend.position = c(0.25, 0.5))

# Add alignment view
n_taxa <- length(alignment)
n_pos <- width(alignment)[1]
trimmedRanges <- GRanges(
  seqnames = names(alignment),
  ranges = IRanges(start = rep(2226, n_taxa), end = rep(n_pos, n_taxa))
)
trimmed_alignment <- alignment[trimmedRanges]
trimmed_alignment <- trimmed_alignment[rotated_tree$tip.label]
bin_alignment <- as.DNAbin(trimmed_alignment)

alignment_tree <- msaplot(p3 + new_scale_fill(), fasta = bin_alignment, 
                          offset = 0.06, width = 0.5) +
  guides(fill = "none") +
  theme(legend.position = c(0.25, 0.5))

# Save output
png_file <- file.path(project_dir, "target_sequencing_tree_alignment.png")
ggsave(alignment_tree, file = png_file, height = 12, width = 8)
cat("Analysis complete. Plot saved as:", png_file, "\n")