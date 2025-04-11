# Swapping Sequence Names in FASTA Alignments with Taxa-Informative Headers

## Overview
This tutorial explains how to replace generic sequence IDs in a FASTA alignment file with more informative taxa names. This is particularly useful for visualizing phylogenetic trees and alignments with meaningful labels.

## Prerequisites
- A FASTA alignment file (`hpc1_aligned.fasta`), from the denovo assembly pipeline
- A mapping file connecting sample numbers to taxa names (`sample_label.tab`), generated using VLOOKUPS from this  [spreadsheet](https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345).

- A mapping file connecting sequence IDs to taxa names (`name_swap.tab`), generated using VLOOKUPS from this  [spreadsheet](https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345).

## Input Files

### 1. Sample Label Mapping (`sample_label.tab`)
This file links the sample S-numbers (from 96-well plates sent for sequencing) to taxa related labels.
See `taxa` field in  this  [spreadsheet](https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345).

```
S60	Bals_S60
S50	Bals_S50
S47	Zlux_S47
S59	Bals_S59
S55	Bals_S55
S62	Bals_S62
S1	Chal_S1
S84	Mesa_S84
S78	Mesa_S78
S14	Chal_S14
```

*Note: This file contains only samples with the introgressed teosinte haplotype for our region of interest.*

### 2. Original FASTA Alignment File
Our starting alignment file (`hpc1_aligned.fasta`) contains headers like:

```
>hpc1_B73  Zm00001eb121780 3:8495640-8499737 dna:chromosome chromosome:Zm-B73-REFERENCE-NAM-5.0:3:1:238017767:1 REF
>hpc1_S29
>hpc1_S82
>hpc1_S26
>hpc1_S41
>hpc1_S72
>hpc1_S45
```

### 3. Name Swap Mapping (`name_swap.tab`)
This file maps sequence IDs to taxa realted names.
See `label` field in  this  [spreadsheet](https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345).

```
hpc1_B73	B73
hpc1_S17	Dura_S17
hpc1_S8	Chal_S8
hpc1_S63	Bals_S63
hpc1_S41	Zlux_S41
hpc1_S83	Bals_S83
hpc1_S10	Chal_S10
hpc1_S50	Bals_S50
hpc1_S30	Zdip_S30
hpc1_S28	Zdip_S28
```

## Step-by-Step Process

### 1. Convert multi-line FASTA to single-line format

```bash
# Convert FASTA with potentially wrapped sequences to one sequence per line
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf "%s",$0}}}' \
  hpc1_aligned.fasta > hpc1_oneline.fasta
```

### 2. Filter the alignment for samples of interest

```bash
# Create empty output file
echo -n > hpc1_filtered.fasta

# Add reference sequences
grep -A1 -p "B73" hpc1_oneline.fasta >> hpc1_filtered.fasta
grep -A1 -p "TIL18" hpc1_oneline.fasta >> hpc1_filtered.fasta

# Add samples from our sample list
cut -f1 sample_label.tab | while IFS=$'\t' read -r SAMPLE; do 
  grep -A1 -p "${SAMPLE}$" hpc1_oneline.fasta >> hpc1_filtered.fasta
done
```

### 3. Re-align the filtered sequences (if needed)

```bash
# Re-align using MAFFT
mafft --reorder --adjustdirection --nuc --auto \
  hpc1_filtered.fasta > hpc1_filtered_realigned.fasta
```

### 4. Replace sequence IDs with taxa names

```bash
# Replace sequence headers using the name_swap.tab mapping file
awk 'FNR==NR{a[$1]=$2;next} /^>/{$0=">"a[substr($0,2)];} 1' \
  name_swap.tab hpc1_filtered_realigned.fasta > hpc1_nice_labels.fasta

# Note: B73 and TIL18 reference sequences may need manual adjustment
```

## Understanding the Name Swapping Command

The key `awk` command breaks down as follows:

1. `FNR==NR{a[$1]=$2;next}` - Read the name mapping file first, storing each pair in the array `a`
2. `/^>/{$0=">"a[substr($0,2)];}` - For each header line (starting with ">"), replace it with ">" followed by the mapped name
3. `substr($0,2)` - Extract the sequence ID by removing the ">" character
4. `1` - Print every line (modified headers and unchanged sequence lines)

## Output Files

After completing these steps, you should have the following files:

```
name_swap/
├── hpc1_aligned.fasta           # Original alignment
├── hpc1_oneline.fasta           # Single-line format of original
├── hpc1_filtered.fasta          # Filtered for samples of interest
├── hpc1_filtered_realigned.fasta # Re-aligned filtered sequences
├── hpc1_nice_labels.fasta       # Final alignment with taxa names
├── name_swap.tab                # Mapping of IDs to taxa names
└── sample_label.tab             # Mapping of sample numbers to labels
```

### Alternative: R Biostrings Approach
For a more readable solution in R, consider using the Biostrings package for FASTA manipulations.
