Maize Genetics: Sequence ID Processing and Labeling Tutorial
================
Your Name
2025-04-14

# Maize Genetics: Sequence ID Processing and Labeling

This tutorial walks through the process of creating standardized
sequence IDs for maize genetic data. We’ll import data from Google
Sheets, process it to create various label formats, and export the
results to a CSV file.

## 1. Introduction

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
project_dir <- "~/Desktop/make_labels"
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

    ## # A tibble: 6 × 16
    ##   genes old_line_id line_id substring sample label accession taxa  `funny_RLLY?`
    ##   <chr> <chr>       <chr>   <chr>     <chr>  <chr> <chr>     <chr>         <dbl>
    ## 1 gpat… JLHNM-661_… Zv.012… Zv0120_P… S60    Bals… Zv.0120   Bals             NA
    ## 2 hpc1… JLNCM-657_… Zv.006… Zv0060_P… S50    Bals… Zv.0060   Bals             NA
    ## 3 hpc1… G6_P4_P4_P… Zl.004… Zl0040_P… S47    Zlux… Zl.0040   Zlux             NA
    ## 4 hpc1… JACV-T-082… Zv.029… Zv0290_P… S59    Bals… Zv.0290   Bals             NA
    ## 5 aaap… JSGYRMM-45… Zv.021… Zv0210_P… S68    Bals… Zv.0210   Bals              1
    ## 6 ipt6… G3_P3_P1_P… Zl.002… Zl0020_P… S40    Zlux… Zl.0020   Zlux             NA
    ## # ℹ 7 more variables: hpc1 <dbl>, gpat15 <dbl>, nrg11 <dbl>, ipt6 <dbl>,
    ## #   tcptf9 <dbl>, nlp15 <dbl>, gdsl <dbl>

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

    ## # A tibble: 6 × 10
    ##   seqid  seqid_sample sample_n maizegdb_prefix taxa  donor_accession donor_sufix
    ##   <chr>  <chr>        <chr>    <chr>           <chr> <chr>           <chr>      
    ## 1 hpc1_… S60          60       Zv              Bals  Zv.0120         0120       
    ## 2 hpc1_… S50          50       Zv              Bals  Zv.0060         0060       
    ## 3 hpc1_… S47          47       Zl              Zlux  Zl.0040         0040       
    ## 4 hpc1_… S59          59       Zv              Bals  Zv.0290         0290       
    ## 5 hpc1_… S68          68       Zv              Bals  Zv.0210         0210       
    ## 6 hpc1_… S40          40       Zl              Zlux  Zl.0020         0020       
    ## # ℹ 3 more variables: ancestry <chr>, ancestry_preffix <chr>, label_1 <chr>

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

    ## # A tibble: 6 × 6
    ##   seqid      founder_ancestry ancestry_preffix label_1  label_2          label_3
    ##   <chr>      <chr>            <chr>            <chr>    <chr>            <chr>  
    ## 1 hpc1_S91   Recurrrent       R                Mesa_S91 R_Zx_Mesa_0510_… Zx_Mes…
    ## 2 hpc1_S34   Recurrrent       R                Nobo_S34 R_Zx_Nobo_0570_… Zx_Nob…
    ## 3 hpc1_S44   Recurrrent       R                Bals_S44 R_Zv_Bals_0440_… Zv_Bal…
    ## 4 hpc1_S49   Recurrrent       R                Bals_S49 R_Zv_Bals_0010_… Zv_Bal…
    ## 5 hpc1_B73   Recurrent        R                B73      R_Zm_B73         Zm_B73 
    ## 6 hpc1_TIL18 Donor            D                TIL18    D_Zx_TIL18       Zx_TIL…

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
    └── sample_annotation.csv
    └── seqid_label.csv

## 6. Data Exploration and Analysis

Let’s explore our final dataset to understand its characteristics:

``` r
# Summary statistics
cat("Total number of entries:", nrow(seqid_label), "\n")
```

    ## Total number of entries: 98

``` r
cat("Distribution by ancestry:\n")
```

    ## Distribution by ancestry:

``` r
print(table(seqid_label$founder_ancestry))
```

    ## 
    ##      Donor  Recurrent Recurrrent 
    ##         53          1         44

``` r
# Example of filtering by ancestry
donor_entries <- seqid_label %>% 
  filter(founder_ancestry == "Donor")

cat("\nDonor entries:", nrow(donor_entries), "\n")
```

    ## 
    ## Donor entries: 53

``` r
knitr::kable(donor_entries, caption = "Donor Entries")
```

| seqid | founder_ancestry | ancestry_preffix | label_1 | label_2 | label_3 |
|:---|:---|:---|:---|:---|:---|
| hpc1_S60 | Donor | D | Bals_S60 | D_Zv_Bals_0120_S60 | Zv_Bals_0120_S60 |
| hpc1_S50 | Donor | D | Bals_S50 | D_Zv_Bals_0060_S50 | Zv_Bals_0060_S50 |
| hpc1_S47 | Donor | D | Zlux_S47 | D_Zl_Zlux_0040_S47 | Zl_Zlux_0040_S47 |
| hpc1_S59 | Donor | D | Bals_S59 | D_Zv_Bals_0290_S59 | Zv_Bals_0290_S59 |
| hpc1_S55 | Donor | D | Bals_S55 | D_Zv_Bals_0420_S55 | Zv_Bals_0420_S55 |
| hpc1_S62 | Donor | D | Bals_S62 | D_Zv_Bals_0260_S62 | Zv_Bals_0260_S62 |
| hpc1_S1 | Donor | D | Chal_S1 | D_Zx_Chal_0160_S1 | Zx_Chal_0160_S1 |
| hpc1_S84 | Donor | D | Mesa_S84 | D_Zx_Mesa_0360_S84 | Zx_Mesa_0360_S84 |
| hpc1_S78 | Donor | D | Mesa_S78 | D_Zx_Mesa_0520_S78 | Zx_Mesa_0520_S78 |
| hpc1_S14 | Donor | D | Chal_S14 | D_Zx_Chal_0030_S14 | Zx_Chal_0030_S14 |
| hpc1_S22 | Donor | D | Dura_S22 | D_Zx_Dura_0550_S22 | Zx_Dura_0550_S22 |
| hpc1_S80 | Donor | D | Mesa_S80 | D_Zx_Mesa_0490_S80 | Zx_Mesa_0490_S80 |
| hpc1_S21 | Donor | D | Dura_S21 | D_Zx_Dura_0560_S21 | Zx_Dura_0560_S21 |
| hpc1_S85 | Donor | D | Bals_S85 | D_Zv_Bals_0150_S85 | Zv_Bals_0150_S85 |
| hpc1_S23 | Donor | D | Mesa_S23 | D_Zx_Mesa_0440_S23 | Zx_Mesa_0440_S23 |
| hpc1_S30 | Donor | D | Zdip_S30 | D_Zd_Zdip_0040_S30 | Zd_Zdip_0040_S30 |
| hpc1_S94 | Donor | D | Mesa_S94 | D_Zx_Mesa_0430_S94 | Zx_Mesa_0430_S94 |
| hpc1_S36 | Donor | D | Zlux_S36 | D_Zl_Zlux_0050_S36 | Zl_Zlux_0050_S36 |
| hpc1_S31 | Donor | D | Nobo_S31 | D_Zx_Nobo_0580_S31 | Zx_Nobo_0580_S31 |
| hpc1_S2 | Donor | D | Chal_S2 | D_Zx_Chal_0160_S2 | Zx_Chal_0160_S2 |
| hpc1_S19 | Donor | D | Dura_S19 | D_Zx_Dura_0550_S19 | Zx_Dura_0550_S19 |
| hpc1_S28 | Donor | D | Zdip_S28 | D_Zd_Zdip_0030_S28 | Zd_Zdip_0030_S28 |
| hpc1_S37 | Donor | D | Hueh_S37 | D_Zh_Hueh_0030_S37 | Zh_Hueh_0030_S37 |
| hpc1_S41 | Donor | D | Zlux_S41 | D_Zl_Zlux_0020_S41 | Zl_Zlux_0020_S41 |
| hpc1_S52 | Donor | D | Bals_S52 | D_Zv_Bals_0340_S52 | Zv_Bals_0340_S52 |
| hpc1_S48 | Donor | D | Bals_S48 | D_Zv_Bals_0370_S48 | Zv_Bals_0370_S48 |
| hpc1_S79 | Donor | D | Bals_S79 | D_Zv_Bals_0480_S79 | Zv_Bals_0480_S79 |
| hpc1_S81 | Donor | D | Mesa_S81 | D_Zx_Mesa_0360_S81 | Zx_Mesa_0360_S81 |
| hpc1_S93 | Donor | D | Mesa_S93 | D_Zx_Mesa_0510_S93 | Zx_Mesa_0510_S93 |
| hpc1_S73 | Donor | D | Zdip_S73 | D_Zd_Zdip_0010_S73 | Zd_Zdip_0010_S73 |
| hpc1_S75 | Donor | D | Bals_S75 | D_Zv_Bals_0310_S75 | Zv_Bals_0310_S75 |
| hpc1_S88 | Donor | D | Mesa_S88 | D_Zx_Mesa_0470_S88 | Zx_Mesa_0470_S88 |
| hpc1_S86 | Donor | D | Mesa_S86 | D_Zx_Mesa_0380_S86 | Zx_Mesa_0380_S86 |
| hpc1_S17 | Donor | D | Dura_S17 | D_Zx_Dura_0560_S17 | Zx_Dura_0560_S17 |
| hpc1_S38 | Donor | D | Hueh_S38 | D_Zh_Hueh_0020_S38 | Zh_Hueh_0020_S38 |
| hpc1_S25 | Donor | D | Dura_S25 | D_Zx_Dura_0540_S25 | Zx_Dura_0540_S25 |
| hpc1_S95 | Donor | D | Mesa_S95 | D_Zx_Mesa_0400_S95 | Zx_Mesa_0400_S95 |
| hpc1_S67 | Donor | D | Zdip_S67 | D_Zd_Zdip_0020_S67 | Zd_Zdip_0020_S67 |
| hpc1_S9 | Donor | D | Chal_S9 | D_Zx_Chal_0140_S9 | Zx_Chal_0140_S9 |
| hpc1_S6 | Donor | D | Chal_S6 | D_Zx_Chal_0190_S6 | Zx_Chal_0190_S6 |
| hpc1_S10 | Donor | D | Chal_S10 | D_Zx_Chal_0340_S10 | Zx_Chal_0340_S10 |
| hpc1_S89 | Donor | D | Mesa_S89 | D_Zx_Mesa_0350_S89 | Zx_Mesa_0350_S89 |
| hpc1_S8 | Donor | D | Chal_S8 | D_Zx_Chal_0090_S8 | Zx_Chal_0090_S8 |
| hpc1_S15 | Donor | D | Mesa_S15 | D_Zx_Mesa_0500_S15 | Zx_Mesa_0500_S15 |
| hpc1_S32 | Donor | D | Nobo_S32 | D_Zx_Nobo_0570_S32 | Zx_Nobo_0570_S32 |
| hpc1_S83 | Donor | D | Bals_S83 | D_Zv_Bals_0470_S83 | Zv_Bals_0470_S83 |
| hpc1_S96 | Donor | D | Mesa_S96 | D_Zx_Mesa_0390_S96 | Zx_Mesa_0390_S96 |
| hpc1_S63 | Donor | D | Bals_S63 | D_Zv_Bals_0140_S63 | Zv_Bals_0140_S63 |
| hpc1_S43 | Donor | D | Bals_S43 | D_Zv_Bals_0030_S43 | Zv_Bals_0030_S43 |
| hpc1_S64 | Donor | D | Bals_S64 | D_Zv_Bals_0220_S64 | Zv_Bals_0220_S64 |
| hpc1_S77 | Donor | D | Bals_S77 | D_Zv_Bals_0250_S77 | Zv_Bals_0250_S77 |
| hpc1_S92 | Donor | D | Mesa_S92 | D_Zx_Mesa_0450_S92 | Zx_Mesa_0450_S92 |
| hpc1_TIL18 | Donor | D | TIL18 | D_Zx_TIL18 | Zx_TIL18 |

Donor Entries

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

    ## R version 4.4.1 (2024-06-14)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS 15.3.2
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] readr_2.1.5 dplyr_1.1.4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.4         knitr_1.50       
    ##  [5] rlang_1.1.6       xfun_0.52         generics_0.1.3    glue_1.8.0       
    ##  [9] bit_4.6.0         htmltools_0.5.8.1 hms_1.1.3         rmarkdown_2.29   
    ## [13] evaluate_1.0.3    tibble_3.2.1      tzdb_0.5.0        fastmap_1.2.0    
    ## [17] yaml_2.3.10       lifecycle_1.0.4   compiler_4.4.1    pkgconfig_2.0.3  
    ## [21] rstudioapi_0.17.1 digest_0.6.37     R6_2.6.1          utf8_1.2.4       
    ## [25] tidyselect_1.2.1  parallel_4.4.1    vroom_1.6.5       pillar_1.10.2    
    ## [29] magrittr_2.0.3    withr_3.0.2       tools_4.4.1       bit64_4.6.0-1
