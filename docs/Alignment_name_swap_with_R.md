# Swapping FASTA Sequence Names in R with Biostrings

This tutorial demonstrates how to use R with the Biostrings package to replace generic sequence IDs in FASTA files with informative taxonomic labels.

## Setup and Required Packages

```r
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
```

## Step 1: Read Mapping Files

First, let's load the mapping tables that connect sequence IDs to taxa names:

```r
# Read the mapping tables
sample_label <- read_tsv("sample_label.tab", col_names = c("sample", "label"))
name_swap <- read_tsv("name_swap.tab", col_names = c("seq_id", "taxa"))

# Display the first few rows of each
head(sample_label)
head(name_swap)
```

## Step 2: Read the FASTA Alignment File

Next, we'll read the original FASTA alignment using Biostrings:

```r
# Read FASTA alignment
original_alignment <- readDNAStringSet("hpc1_aligned.fasta")

# View some information about the alignment
length(original_alignment)  # Number of sequences
width(original_alignment)   # Sequence lengths
names(original_alignment)[1:10]  # First 10 sequence names
```

## Step 3: Filter Sequences of Interest

Let's filter the alignment to keep only our samples of interest:

```r
# Create a vector of sequences to keep (reference sequences + samples)
samples_to_keep <- c("hpc1_B73", "hpc1_TIL18", paste0("hpc1_", sample_label$sample))

# Filter the alignment
filtered_alignment <- original_alignment[names(original_alignment) %in% samples_to_keep]

# Check how many sequences we have now
length(filtered_alignment)
```

## Step 4: Create a Complete Name Mapping Table

Now we'll build a comprehensive mapping table for all sequence IDs:

```r
# Extract sequence IDs from the filtered alignment
seq_ids <- names(filtered_alignment)

# Strip "hpc1_" prefix for matching
simple_ids <- gsub("^hpc1_", "", seq_ids)

# Create a comprehensive mapping dataframe
name_mapping <- data.frame(
  seq_id = seq_ids,
  stringsAsFactors = FALSE
)

# Add taxa names using the mapping tables
name_mapping <- name_mapping %>%
  left_join(name_swap, by = "seq_id") %>%
  # For any missing mappings, try to match with sample_label
  mutate(
    simple_id = gsub("^hpc1_", "", seq_id),
    taxa = ifelse(is.na(taxa) & simple_id %in% sample_label$sample,
                  sample_label$label[match(simple_id, sample_label$sample)],
                  taxa)
  ) %>%
  # For any remaining NAs, keep the original name
  mutate(taxa = ifelse(is.na(taxa), seq_id, taxa))

# View the mapping table
head(name_mapping)
```

## Step 5: Replace Sequence Names and Write New FASTA

Finally, we'll apply the mapping and write out the new alignment:

```r
# Create a named vector for easy lookup
name_map_vector <- setNames(name_mapping$taxa, name_mapping$seq_id)

# Apply the name changes
names(filtered_alignment) <- name_map_vector[names(filtered_alignment)]

# Write the renamed alignment to a new file
writeXStringSet(filtered_alignment, "hpc1_nice_labels.fasta", format="fasta")
```

## Step 6 (Optional): Realign Sequences with R

If you need to re-align your sequences, you can use the DECIPHER package:

```r
# Install and load DECIPHER if needed
if (!requireNamespace("DECIPHER", quietly = TRUE))
  BiocManager::install("DECIPHER")
library(DECIPHER)

# Perform multiple sequence alignment
aligned_sequences <- AlignSeqs(filtered_alignment, iterations = 2)

# Write the aligned sequences with nice labels
writeXStringSet(aligned_sequences, "hpc1_realigned_nice_labels.fasta", format="fasta")
```

## Complete Script

Here's the complete R script combining all steps:

```r
# Load required packages
library(Biostrings)
library(tidyverse)

# Read mapping tables
sample_label <- read_tsv("sample_label.tab", col_names = c("sample", "label"))
name_swap <- read_tsv("name_swap.tab", col_names = c("seq_id", "taxa"))

# Read the original alignment
original_alignment <- readDNAStringSet("hpc1_aligned.fasta")

# Create a vector of sequences to keep
samples_to_keep <- c("hpc1_B73", "hpc1_TIL18", paste0("hpc1_", sample_label$sample))

# Filter the alignment
filtered_alignment <- original_alignment[names(original_alignment) %in% samples_to_keep]

# Create a comprehensive mapping dataframe
name_mapping <- data.frame(
  seq_id = names(filtered_alignment),
  stringsAsFactors = FALSE
) %>%
  left_join(name_swap, by = "seq_id") %>%
  mutate(
    simple_id = gsub("^hpc1_", "", seq_id),
    taxa = ifelse(is.na(taxa) & simple_id %in% sample_label$sample,
                  sample_label$label[match(simple_id, sample_label$sample)],
                  taxa)
  ) %>%
  mutate(taxa = ifelse(is.na(taxa), seq_id, taxa))

# Create a named vector for easy lookup
name_map_vector <- setNames(name_mapping$taxa, name_mapping$seq_id)

# Apply the name changes
names(filtered_alignment) <- name_map_vector[names(filtered_alignment)]

# Write the renamed alignment to a new file
writeXStringSet(filtered_alignment, "hpc1_nice_labels.fasta", format="fasta")

# Optional: Realign with DECIPHER
# if (requireNamespace("DECIPHER", quietly = TRUE)) {
#   library(DECIPHER)
#   aligned_sequences <- AlignSeqs(filtered_alignment, iterations = 2)
#   writeXStringSet(aligned_sequences, "hpc1_realigned_nice_labels.fasta", format="fasta")
# }
```

## Advantages of the R Approach

1. **Readability**: The R code is more self-explanatory than complex bash commands
2. **Flexibility**: Easy to modify for different naming patterns or additional data processing
3. **Error Handling**: Better handling of edge cases and missing mappings
4. **Integration**: Can be part of a larger R-based bioinformatics workflow
5. **Visualization**: Can easily add sequence visualization or analysis in the same script

This approach eliminates the need for intermediate files and complex awk or grep commands, making the process more streamlined and reproducible.