# Required libraries
library(targseq)
library(Biostrings)
library(phangorn)
library(ape)
library(ggtree)
library(ggplot2)
library(colorspace)
library(dplyr)

# 1. Read data
aln_file <- system.file("extdata", "hpc1_aligned.fasta", package="targseq")
alignment <- readDNAStringSet(aln_file)
alignment_order <- data.frame(seqid=names(alignment)) 

# metadata_file <- system.file("extdata", "seqid_label.csv", package="targseq")
metadata_file <- ('~/Desktop/seqid_label.csv')
donor_file <- ('~/Desktop/donor_data.csv')
donor_data <- read.csv(donor_file)

# Read metadata 
cat("Reading sample metadata...\n")

taxa_info <- read.csv(metadata_file)

taxa_info <- alignment_order  %>%
  left_join(taxa_info) %>%           
  left_join(donor_data, by=c(donor_accession= "donor_id")) %>%
  select(seqid, label=label_3,ancestry_call,elevation, nitrogen_0.5cm_mean_1000)


# taxa_info <- read.csv(system.file("extdata", "seqid_label.csv", package="targseq"))
rownames(taxa_info) <- NULL


# 2. Extract variant genotypes
variants <- data.frame(
  name = c("A204T", "I211V"),
  pattern = c("GCCGTGGCGTGGCGC", "ATCACCCGC")
)
variant_data <- get_variant_gt(variants, alignment)

# 3. Recode variants as REF/ALT
variant_data$A204T <- c("ALT", "REF")[as.factor(variant_data$A204T)]
variant_data$I211V <- c("REF", "ALT")[as.factor(variant_data$I211V)]
rownames(variant_data) <-taxa_info$label

# 4. Build tree
phyDat <- phyDat(read.dna(aln_file, format="fasta"), type="DNA")
names(phyDat) <- taxa_info$label
dna_dist <- dist.ml(phyDat, model="JC69")
tree <- ladderize(upgma(dna_dist), right=FALSE)
rotated_tree <- pivot_on(tree, "Zm_B73")

# 5. Create visualization
tree_plot <- create_variant_heatmap_tree(
  tree = rotated_tree, 
  data = variant_data
)
taxa_info$ancestry_call
taxa_info$label
rotated_tree$tip.label
hist(taxa_info$nitrogen_0.5cm_mean_1000)

ggtree(rotated_tree,ladderize = FALSE) %<+% taxa_info[,-1] +
  geom_tippoint(aes(color=nitrogen_0.5cm_mean_1000)) +
  scale_color_continuous_divergingx(palette = 'ArmyRose', mid =300, min=100,n_interp = 25)  +
  theme(legend.position = c(0.25, 0.5))

ggtree(rotated_tree,ladderize = FALSE) %<+% taxa_info[,-1] +
    geom_tippoint(aes(color=elevation)) +
  scale_color_continuous_divergingx(palette = 'RdBu', mid =1500, n_interp = 25)  +
  theme(legend.position = c(0.25, 0.5))
                


# 6. Add ancestry information
ancestry_palette <- c("Recurrent" = "tomato", "Donor" = "royalblue")
enhanced_tree <- tree_plot %<+% taxa_info[,-1] +
  geom_tippoint(aes(color = ancestry_call), position = position_nudge(x = 0.0015)) +
   scale_color_manual(values = ancestry_palette) +
  theme(legend.position = c(0.25, 0.5))


# 7. Add alignment view
alignment_tree <- msaplot(enhanced_tree, 
                          fasta = aln_file, 
                          offset = 0.05, 
                          width = 0.6)

# 8. Save outputs
ggsave(enhanced_tree, file = "maize_tree.pdf", height = 10, width = 7)
ggsave(alignment_tree, file = "maize_tree_alignment.png", height = 12, width = 8)
