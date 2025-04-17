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

line_info_file <- file.path(project_dir, "FINAL_SELECTION.csv")
sample_annotation_file <- file.path(project_dir, "sample_annotation.csv")

# Import data from local CSV files
line_info <- read_csv(line_info_file)
sample_annotation <- read_csv(sample_annotation_file)

# Display the first few rows of line_info
print("Input data preview:")
print(head(line_info))

# Create sequence ID table with informative labels
seqid_label <- line_info %>%
  # Select and rename columns for clarity
  select(
    sample,                     # Sample identifier
    donor_accession = accession, # Rename for clarity
    taxa,                       # Taxonomic information
    introgression = hpc1,       # Introgression status (0/1)
    label_1 = label             # Original label
  ) %>%
  # Create new metadata columns through transformations
  mutate(
    # Extract prefix from donor accession (before the first period)
    maizegdb_prefix = gsub("\\..*", "", donor_accession, perl = TRUE),
    
    # Create standardized sequence ID by combining "hpc1" and sample
    seqid = paste("hpc1", sample, sep = "_"),
    
    # Copy sample ID for reference
    seqid_sample = sample,
    
    # Remove "S" prefix from sample numbers
    sample_n = gsub("S", "", sample, fixed = TRUE),
    
    # Extract suffix from donor accession (after the last period)
    donor_sufix = gsub(".*\\.", "", donor_accession, perl = TRUE),
    
    # Create ancestry prefix: "R" for Recurrent (introgression=0), "D" for Donor (introgression=1)
    ancestry_preffix = c("R", "D")[introgression + 1],
    
    # Create full ancestry label: "Recurrent" or "Donor"
    ancestry = c("Recurrrent", "Donor")[introgression + 1]
  ) %>%
  # Select and reorder columns for final output
  select(
    seqid, seqid_sample, sample_n, maizegdb_prefix, taxa, 
    donor_accession, taxa, donor_sufix, ancestry, ancestry_preffix, label_1
  )

# Display the initial processing results
print("Initial sequence ID table:")
print(head(seqid_label))

# Add additional label formats and reference genomes
seqid_label <- rbind(
  # Process existing entries with additional labels
  seqid_label %>%
    mutate(
      # Create label_2: Combined format with ancestry prefix
      # This keeps taxonomic order from general to specific
      label_2 = paste(ancestry_preffix, maizegdb_prefix, taxa, donor_sufix, seqid_sample, sep = "_"),
      
      # Create label_3: Format without ancestry prefix
      label_3 = paste(maizegdb_prefix, taxa, donor_sufix, seqid_sample, sep = "_")
    ) %>% 
    select(seqid, founder_ancestry = ancestry, ancestry_preffix, starts_with("label")),
  
  # Add B73 reference genome entry (Recurrent)
  data.frame(
    seqid = "hpc1_B73", 
    founder_ancestry = "Recurrent", 
    ancestry_preffix = "R",
    label_1 = "B73", 
    label_2 = "R_Zm_B73", 
    label_3 = "Zm_B73",
    stringsAsFactors = FALSE
  ),
  
  # Add TIL18 reference genome entry (Donor)
  data.frame(
    seqid = "hpc1_TIL18", 
    founder_ancestry = "Donor", 
    ancestry_preffix = "D",
    label_1 = "TIL18", 
    label_2 = "D_Zx_TIL18", 
    label_3 = "Zx_TIL18",
    stringsAsFactors = FALSE
  )
)

# Display the final table
print("Final sequence ID table with reference genomes:")
print(tail(seqid_label))

# Set output file path
output_file <- file.path(project_dir, "seqid_label_table.csv")

# Write the table to CSV
write.csv(seqid_label, output_file, quote = FALSE, row.names = FALSE)

cat("Sequence ID table successfully created and saved to:", output_file, "\n")

# Data Exploration and Analysis
cat("\n--- Data Exploration ---\n")
# Summary statistics
cat("Total number of entries:", nrow(seqid_label), "\n")
cat("Distribution by ancestry:\n")
print(table(seqid_label$founder_ancestry))

# Example of filtering by ancestry
donor_entries <- seqid_label %>% 
  filter(founder_ancestry == "Donor")

cat("\nDonor entries:", nrow(donor_entries), "\n")
print(head(donor_entries))

cat("\nScript execution complete.\n")
