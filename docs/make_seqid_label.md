Maize Genetics: Sequence ID Processing and Labeling Tutorial
================
Your Name
2025-04-14

## 1. Introduction

This tutorial walks through the process of creating standardized
sequence IDs for maize genetic data. We’ll import data from Google
Sheets, process it to create various label formats, and export the
results to a CSV file.

Consistent labeling of sequence IDs is crucial for genetic data
analysis, especially when working with introgression lines or
comparative genomics. This script demonstrates how to create a
standardized naming system that encodes important metadata directly in
the sequence IDs.

### Starting Directory Structure

Before running the script, your project directory should look like this:

    make_labels/
    ├── FINAL_SELECTION.csv
    └── sample_annotation.csv

``` r
# Load required libraries
library(dplyr)    # For data manipulation
library(readr)    # For reading/writing CSV files
```

Now we’ll import the data from local CSV files. Make sure to download
these files from Google Sheets first and save them to your working
directory.

``` r
# Create a project directory if it doesn't exist
project_dir <- "make_labels"
if (!dir.exists(project_dir)) {
  dir.create(project_dir)
}

# Set file paths within the project directory
line_info_file <- file.path(project_dir, "FINAL_SELECTION.csv")
sample_annotation_file <- file.path(project_dir, "sample_annotation.csv")


# Import data from local CSV files
# Make sure these files are in your working directory
line_info <- read_csv(line_info_file)
sample_annotation <- read_csv(sample_annotation_file)

# Display the first few rows of line_info
head(line_info)
```

## 2. Create sequence ID table with informative labels

This step selects relevant columns, extracts components from accession
IDs, and creates standardized sequence IDs with useful metadata

``` r
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
head(seqid_label)
```

Let’s examine what happened in the code above:

1.  We selected specific columns from the input data
2.  We extracted the prefix from accessions (e.g., “Zm” from “Zm.B73.1”)
3.  We created standardized sequence IDs by prefixing “hpc1\_” to sample
    IDs
4.  We extracted numeric sample values (removing the “S” prefix)
5.  We determined ancestry categories (Recurrent/Donor) based on
    introgression status

## 4. Adding Additional Label Formats and Reference Genomes

Now we’ll add additional label formats and include reference genome
entries:

``` r
# Add additional label formats and reference genomes
seqid_label <- rbind(
  # Process existing entries with additional labels
  seqid_label %>%
    mutate(
      # Create label_2: Combined format with ancestry prefix
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
tail(seqid_label)
```

In this section:

1.  We created two additional label formats:
    - `label_2`: Includes ancestry prefix, MaizeGDB prefix, taxa,
      suffix, and sample ID. ***This is in order of taxonomic ranks, it
      keeps taxonomic order from general to specific***
    - `label_3`: Similar to label_2 but without the ancestry prefix
2.  We added entries for two reference genomes:
    - B73 (Recurrent parent)
    - TIL18 (Donor parent)

## 5. Exporting the Results

Now we’ll export our processed data to a CSV file:

``` r
output_file <- file.path(project_dir, "seqid_label_table.csv")

# Write the table to CSV
write.csv(seqid_label, output_file, quote = FALSE, row.names = FALSE)

cat("Sequence ID table successfully created and saved to:", output_file, "\n")
```

    make_labels/
    ├── FINAL_SELECTION.csv
    ├── sample_annotation.csv
    └── seqid_label.csv

## 6. Data Exploration and Analysis

Let’s explore our final dataset to understand its characteristics:

``` r
# Summary statistics
cat("Total number of entries:", nrow(seqid_label), "\n")
cat("Distribution by ancestry:\n")
print(table(seqid_label$founder_ancestry))

# Example of filtering by ancestry
donor_entries <- seqid_label %>% 
  filter(founder_ancestry == "Donor")

cat("\nDonor entries:", nrow(donor_entries), "\n")
head(donor_entries)
```

## 7. Conclusion

In this tutorial, we’ve demonstrated how to:

1.  Import maize genetic data (from Google Sheets or local files)
2.  Process the data to create standardized sequence IDs
3.  Generate multiple label formats for different purposes
4.  Add reference genome entries
5.  Export the results to a CSV file
6.  Perform basic data exploration

This standardized labeling system makes it easier to organize and
analyze genetic data, especially when working with introgression lines
or comparative genomics projects.

## 8. Session Information

For reproducibility, here’s the session information:

``` r
sessionInfo()
```
