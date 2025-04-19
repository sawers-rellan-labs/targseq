# Required libraries
library(targseq)
library(Biostrings)
library(phangorn)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

# 1. Read data
aln_file <- system.file("extdata", "hpc1_aligned.fasta", package="targseq")
alignment <- readDNAStringSet(aln_file)
taxa_info <- read.csv(system.file("extdata", "seqid_label.csv", package="targseq"))
rownames(taxa_info) <- taxa_info$seqid
taxa_info <- taxa_info %>% rename(ancestry_call = "founder_ancestry")

# 2. Extract variant genotypes
variants <- data.frame(
  name = c("A204T", "I211V"),
  pattern = c("GCCGTGGCGTGGCGC", "ATCACCCGC")
)
variant_data <- get_variant_gt(variants, alignment)
variant_data$A204T <- c("ALT", "REF")[as.factor(variant_data$A204T)]
variant_data$I211V <- c("REF", "ALT")[as.factor(variant_data$I211V)]
rownames(variant_data) <- names(alignment)

# 3. Build tree
phyDat <- phyDat(read.dna(aln_file, format="fasta"), type="DNA")
dna_dist <- dist.ml(phyDat, model="JC69")
tree <- ladderize(upgma(dna_dist), right=FALSE)
rotated_tree <- pivot_on(tree, "hpc1_B73")

# 4. Create visualization
tree_plot <- create_variant_heatmap_tree(
  tree = rotated_tree, 
  data = variant_data
)

# 5. Add ancestry information
ancestry_palette <- c("Recurrent" = "tomato", "Donor" = "royalblue")
enhanced_tree <- tree_plot %<+% taxa_info +
  geom_tippoint(aes(color = ancestry_call), position = position_nudge(x = 0.0015)) +
  scale_color_manual(values = ancestry_palette) +
  theme(legend.position = c(0.25, 0.5))

# 6. Add alignment view
alignment_tree <- msaplot(enhanced_tree, 
                          fasta = aln_file, 
                          offset = 0.05, 
                          width = 0.6)

# 7. Save outputs
ggsave(enhanced_tree, file = "maize_tree.pdf", height = 10, width = 7)
ggsave(alignment_tree, file = "maize_tree_alignment.png", height = 12, width = 8)
