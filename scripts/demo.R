# Installation steps for the targseq package
# Run these steps first to install the package from source

# Install required development tools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# install from GitHub if you've pushed the package there
devtools::install_github("sawers-rellan-labs/targseq")

# Example usage of the targseq package
library(targseq)
library(Biostrings)
library(phangorn)
library(ape)
library(ggtree)
library(ggplot2)

# Example workflow
# 1. Read alignment file

aln_file <- system.file("extdata", "hpc1_aligned.fasta", package="targseq")

alignment <- readDNAStringSet(aln_file)

# 2. Extract variant information
variants <- data.frame(
  name = c("A204T", "I211V"),
  pattern = c("GCCGTGGCGTGGCGC", "ATCACCCGC")
)

variant_data <- get_variant_gt(variants,alignment)

# 3. Build phylogenetic tree
phyDat  <- phyDat(data = read.dna(aln_file, format="fasta"), type = "DNA", levels = NULL)

# Calculate distance matrix using JC69 model
# Note: You could use modelTest() to find the best model for your data
# mt <- modelTest(hpc1_phyDat)
dna_dist <- dist.ml(phyDat, model="JC69")
tree <- ladderize(upgma(dna_dist), right = FALSE)
# 4. Rotate tree to put reference at top
rotated_tree <- pivot_on(tree, "hpc1_B73")

# 5. Visualize tree with variants

# Change coding of the genotype to REF/ALT  
variant_data$A204T<- c("ALT","REF")[as.factor(variant_data$A204T)]
variant_data$I211V<- c("REF","ALT")[as.factor(variant_data$I211V)]
rownames(variant_data)<- names(alignment)


taxa_info <- read.csv(system.file("extdata", "seqid_label.csv", package="targseq"))
rownames(taxa_info) <-taxa_info$seqid
taxa_info <- taxa_info[tree$tip.label,] %>%
  rename(founder_ancestry="ancestry_call")
# Create visualization

tree_plot <- create_variant_heatmap_tree(
  tree = rotated_tree, 
  data = variant_data)

quartz(height=12); print(tree_plot)

pal <-c("Recurrent" = "tomato", "Donor" = "royalblue")

quartz()
tree_plot <- tree_plot %<+% taxa_info  +
  geom_tippoint(aes(color = ancestry_call ), position = position_nudge(x = 0.0015)) +
  scale_color_manual(values=pal) +
  theme(legend.position = c(0.25,0.5))


quartz(height=12); print(tree_plot)

p <- ggtree(rotated_tree, ladderize = FALSE) %<+% 
  taxa_info  + 
  geom_tiplab(size = 3) +
  geom_tippoint(aes(color = ancestry_call ), position = position_nudge(x = 0.0015)) +
  scale_color_manual(values=pal)

p2 <- msaplot(tree_plot, fasta = aln_file, 
              offset = 0.05, width = 0.6)  +
  theme(legend.position = c(0.25,0.5))


quartz(height=12); print(p2)

# Save the tree with alignment visualization
ggsave(p2, file = file.path(project_dir, "hpc1_tree_alignment.png"), 
       height = 12, width = 7, units = "in")

