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


taxa_info <- system.file("extdata", "seqid_label.csv", package="targseq")

# Create visualization
tree_plot <- create_variant_heatmap_tree(
  tree = rotated_tree, 
  data =variant_data)

tree_plot %<+% variant_data +
  geom_tippoint(aes(color=Country)) 

print(tree_plot)
