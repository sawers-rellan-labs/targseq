#!/usr/bin/env Rscript

# Maize Genetics: Sequence ID Processing and Labeling
# Script to create standardized sequence IDs for maize genetic data
# Extracted from R Markdown tutorial

# Load required libraries
library(dplyr)    # For data manipulation
library(readr)    # For reading/writing CSV files

# Create a project directory if it doesn't exist
project_dir <- "make_labels"
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}

# Set file paths within the project directory
# data from
# https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345
# Load necessary libraries
library(dplyr)
library(tidyr)

# Read input files
final_selection <- read.csv("~/Desktop/FINAL_SELECTION.csv")

sample_annotation <- read.csv("~/Desktop/sample_annotation.csv") %>%
  select(fastq_prefix,fastq_sample) %>%
  arrange(fastq_prefix,fastq_sample) %>%
  distinct() %>%
  # Remove "S" prefix from sample numbers
  mutate(sample_n = gsub("S","",fastq_sample, fixed=TRUE) %>% as.numeric()) %>%
  arrange(sample_n)

target_genes <- c("hpc1", "gpat15", "nrg11", "ipt6", "tcptf9", "nlp15", "gdsl")
# Join the datasets
seqid_label_wide <- final_selection%>%
  left_join(sample_annotation, by = c("fastq_prefix", "fastq_sample"))

# Pivot the data from wide to long format
seqid_label <- seqid_label_wide  %>%
  dplyr::select(-genes, -hpc1_ancestry_mismatch) %>%
  pivot_longer(
    cols = all_of(target_genes),
    names_to = "gene",
    values_to = "ancestry_call"
  ) %>%
  select(gene,sample_n,fastq_sample, everything()) %>%
  arrange(gene,sample_n)

colnames(seqid_label)

seqid_label <- seqid_label %>%
  # Create new metadata columns through transformations
  mutate(
    # Extract prefix from donor accession (before the first period)
    maizegdb_prefix = gsub("\\..*", "", donor_accession, perl = TRUE),
    
    # Create standardized sequence ID by combining "hpc1" and sample
    seqid = paste(gene, fastq_sample, sep = "_"),
    
    label_1 = seqid,
    
    # Copy sample ID for reference
    seqid_sample = fastq_sample,
    
    # Extract suffix from donor accession (after the last period)
    donor_suffix = gsub(".*\\.", "", donor_accession, perl = TRUE),
    
    ancestry_call =  factor(ancestry_call,levels= c("Recurrent","Donor")),
    
    # Create ancestry prefix: "R" for Recurrent (introgression=0), "D" for Donor (introgression=1)
    ancestry_prefix = c("R", "D")[ancestry_call],
  ) %>%
  # Select and reorder columns for final output
  select(
    gene,seqid, seqid_sample, sample_n, maizegdb_prefix, taxa, 
    donor_accession, taxa, donor_suffix, ancestry_call, ancestry_prefix, label_1
  )



ref_string <- "gene sample_n seqid  ancestry_call  ancestry_prefix donor_accession label_1 label_2 label_3
              hpc1 NA hpc1_B73  Recurrent  R NA  hpc1_B73 R_Zm_B73 Zm_B73
              hpc1 NA  hpc1_TIL18  Donor  D NA hpc1_TIL18 D_Zx_TIL18 Zx_TIL18"

ref_label <- lapply(target_genes, FUN = function(x){
  data_string <- gsub("hpc1", x,ref_string, fixed= TRUE) 
  read.table(text =data_string, header = TRUE, na.strings = "NA")
}) %>%  bind_rows()


# Display the initial processing results
head(seqid_label)


# Add additional label formats and reference genomes
seqid_label <-  rbind ( seqid_label %>%
                          # Process existing entries with additional labels
                          
                          mutate(
                            # Create label_2: Combined format with ancestry prefix
                            label_1 = seqid,
                            
                            # Create label_2: Combined format with ancestry prefix
                            label_2 = paste(ancestry_prefix, maizegdb_prefix, taxa, donor_suffix, seqid_sample, sep = "_"),
                            
                            # Create label_3: Format without ancestry prefix
                            label_3 = paste(maizegdb_prefix, taxa, donor_suffix, seqid_sample, sep = "_")
                          ) %>% 
                          select(gene,sample_n,seqid, ancestry_call, ancestry_prefix, donor_accession, starts_with("label")) %>%
                          group_by(gene) %>%
                          arrange(gene,sample_n),
                        ref_label)

write.csv(seqid_label,file="~/Desktop/seqid_label.csv",row.names = FALSE)

# Display the final table
tail(seqid_label)

