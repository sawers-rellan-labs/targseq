#!/usr/bin/env Rscript

#===============================================================================
# Maize Genetics: Sequence ID Processing and Labeling
# Author: [Your Name]
# Date: [Current Date]
#
# Description:
#   This script processes maize genetic data to create standardized sequence IDs.
#   It combines data from the FINAL_SELECTION.csv with sample annotations to 
#   generate formatted sequence IDs with ancestry information for target genes.
#
# Input Files:
#   - FINAL_SELECTION.csv: Contains gene line IDs and ancestry calls
#   - sample_annotation.csv: Maps fastq prefixes to sample IDs
#
# Output:
#   - seqid_label.csv: Long-format table with standardized sequence IDs and labels
#===============================================================================

#-------------------------------------------------------------------------------
# 1. SETUP AND CONFIGURATION
#-------------------------------------------------------------------------------

# Load required libraries
library(dplyr)    # For data manipulation
library(readr)    # For reading/writing CSV files
library(tidyr)    # For data reshaping

# Create a project directory if it doesn't exist
project_dir <- "make_labels"
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}

# Define genes for target sequencing  
target_genes <- c("hpc1", "gpat15", "nrg11", "ipt6", "tcptf9", "nlp15", "gdsl")

#-------------------------------------------------------------------------------
# 2. DATA IMPORT
#-------------------------------------------------------------------------------

# Read the main selection data
# Source: https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit

final_selection <- read.csv("~/Desktop/FINAL_SELECTION.csv")

# Read and preprocess sample annotation data

sample_annotation <- read.csv("~/Desktop/sample_annotation.csv") %>%
  select(fastq_prefix, fastq_sample) %>%
  arrange(fastq_prefix, fastq_sample) %>%
  distinct() %>%
  # Remove "S" prefix from sample numbers and convert to numeric for sorting
  mutate(sample_n = gsub("S", "", fastq_sample, fixed=TRUE) %>% as.numeric()) %>%
  arrange(sample_n)

#-------------------------------------------------------------------------------
# 3. DATA JOINING AND TRANSFORMATION
#-------------------------------------------------------------------------------

# Join datasets to combine sample information with genetic data
seqid_label_wide <- final_selection %>%
  left_join(sample_annotation)

# Transform data from wide to long format (one row per gene-sample combination)
seqid_label <- seqid_label_wide  %>%
  # Remove unnecessary columns
  dplyr::select(-genes) %>%
  # Pivot gene columns to long format
  pivot_longer(
    cols = all_of(target_genes),
    names_to = "gene",
    values_to = "ancestry_call"
  ) %>%
  # Reorganize columns and sort by gene and sample number
  select(gene, sample_n, fastq_sample, everything()) %>%
  arrange(gene, sample_n)

# Display column names for verification
colnames(seqid_label)

#-------------------------------------------------------------------------------
# 4. SEQUENCE ID AND LABEL GENERATION
#-------------------------------------------------------------------------------

# Create standardized labels and metadata
seqid_label <- seqid_label %>%
  # Generate additional metadata columns
  mutate(
    # Extract prefix from donor accession (before the first period)
    maizegdb_prefix = gsub("\\..*", "", donor_accession, perl = TRUE),
    
    # Create standardized sequence ID by combining gene and sample
    seqid = paste(gene, fastq_sample, sep = "_"),
    
    # Initialize primary label
    label_1 = seqid,
    
    # Store sample ID for reference
    seqid_sample = fastq_sample,
    
    # Extract suffix from donor accession (after the last period)
    donor_suffix = gsub(".*\\.", "", donor_accession, perl = TRUE),
    
    # Convert ancestry call to factor with ordered levels
    ancestry_call = factor(ancestry_call, levels = c("Recurrent", "Donor")),
    
    # Create ancestry prefix: "R" for Recurrent, "D" for Donor
    ancestry_prefix = c("R", "D")[ancestry_call],
  ) %>%
  # Select and reorder columns for intermediate output
  select(
    gene, fastq_prefix, fastq_sample, seqid, seqid_sample, sample_n, maizegdb_prefix, taxa, 
    donor_accession, taxa, donor_suffix, ancestry_call, ancestry_prefix, label_1
  )

# Display the intermediate processing results
head(seqid_label)

#-------------------------------------------------------------------------------
# 5. REFERENCE DATA CREATION
#-------------------------------------------------------------------------------

# Taxa database table
taxa_db <- read.table(text = "
B73     Zm-B73-REFERENCE-NAM-5.0   1
TIL18   Zx-TIL18-REFERENCE-PanAnd-1.0      2
TIL01   Zv-TIL01-REFERENCE-PanAnd-1.0      3
TdFL    Td-FL_9056069_6-REFERENCE-PanAnd-2.0a      4
Zd      Zd-Gigi-REFERENCE-PanAnd-1.0       5
Zh      Zh-RIMHU001-REFERENCE-PanAnd-1.0   6
PT      Zm-PT-REFERENCE-HiLo-1.0   7
TIL11   Zv-TIL11-REFERENCE-PanAnd-1.0      8
TIL25   Zx-TIL25-REFERENCE-PanAnd-1.0      9
Momo    Zd-Momo-REFERENCE-PanAnd-1.0       10
Zn      Zn-PI615697-REFERENCE-PanAnd-1.0   11", 
                      header = FALSE, 
                      col.names = c("accession", "reference_name", "id"))

# Extract species prefix from reference_name
taxa_db <- taxa_db %>%
  mutate(species_prefix = sub("^(\\w+)-.*$", "\\1", reference_name),
         accession_name = sub("^\\w+-([^-]+).*$", "\\1", reference_name))

# List of target genes
target_genes <- c("ipt6", "hpc1", "nrg11", "nlp1", "gdsl", "tcptf9")
gene_ids <- c(
  "ipt6" = "Zm00001eb062030",
  "hpc1" = "Zm00001eb121780",
  "nrg11" = "Zm00001eb206940",
  "nlp1" = "Zm00001eb231720",
  "gdsl" = "Zm00001eb268440",
  "tcptf9" = "Zm00001eb372490"
)

# Create a data frame with all combinations of genes and accessions
ref_label<- expand.grid(
  gene = target_genes,
  accession = taxa_db$accession,
  stringsAsFactors = FALSE
) %>%
  left_join(taxa_db, by = "accession") %>%
  mutate(
    fastq_prefix = NA,
    fastq_sample = NA,
    sample_n = NA,
    seqid = paste0(gene, "_", accession),
    ancestry_call = case_when(
      accession == "B73" ~ "Recurrent",
      accession == "TdFl" ~ NA,
      .default ="Donor"),
    ancestry_prefix = case_when(
      accession == "B73" ~ "R", 
      accession == "TdFL" ~ NA,
      .default = "D"),
    donor_accession = NA,
    label_1 = seqid
  )

# Fix the redundancy issue in label_2 and label_3 for Zd, Zn, and Zh
ref_label <- ref_label %>%
  mutate(
    # Handle the special cases for label_2
    label_2 = case_when(
      accession == "Zd" & species_prefix == "Zd" ~ paste0(ancestry_prefix, "_Zd"),
      accession == "Zn" & species_prefix == "Zn" ~ paste0(ancestry_prefix, "_Zn"),
      accession == "Zh" & species_prefix == "Zh" ~ paste0(ancestry_prefix, "_Zh"),
      TRUE ~ paste0(ancestry_prefix, "_", species_prefix, "_", accession)
    ),
    # Handle the special cases for label_3
    label_3 = case_when(
      accession == "Zd" & species_prefix == "Zd" ~ "Zd",
      accession == "Zn" & species_prefix == "Zn" ~ "Zn",
      accession == "Zh" & species_prefix == "Zh" ~ "Zh",
      TRUE ~ paste0(species_prefix, "_", accession)
    )
  ) %>%
  select(
    gene, fastq_prefix, fastq_sample, sample_n, seqid, 
    ancestry_call, ancestry_prefix, donor_accession, 
    label_1, label_2, label_3
  )

# View the result
head(ref_label)
ref_label

# Write the result to a file if needed
# write.table(result, "updated_gene_table.tsv", sep="\t", quote=FALSE, row.names=FALSE)


#-------------------------------------------------------------------------------
# 6. FINAL DATA ASSEMBLY AND EXPORT
#-------------------------------------------------------------------------------

# Combine processed data with reference data and add additional label formats
seqid_label <- rbind(
  seqid_label %>%
    # Generate additional label formats
    mutate(
      # Ensure basic label is set
      label_1 = seqid,
      
      # Create label_2: Format with ancestry prefix
      label_2 = paste(ancestry_prefix, maizegdb_prefix, taxa, donor_suffix, seqid_sample, sep = "_"),
      
      # Create label_3: Format without ancestry prefix
      label_3 = paste(maizegdb_prefix, taxa, donor_suffix, seqid_sample, sep = "_")
    ) %>% 
    # Select final columns and sort
    select(gene, fastq_prefix,fastq_sample, sample_n, seqid, ancestry_call, ancestry_prefix, donor_accession, starts_with("label")) %>%
    group_by(gene) %>%
    arrange(gene, sample_n),
  
  # Add reference data
  ref_label
)

# Write the final data to CSV
write.csv(seqid_label, file = "~/Desktop/seqid_label.csv", quote=FALSE, na='',row.names = FALSE)

# Display the final table
head(seqid_label)
tail(seqid_label)
