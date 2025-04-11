# Swapping sequence names in  the fasta alignment with a taxa informative header

The `sample_label.tab` links the sample `S` number, from the 96 well plate sent to sequencing, to the `label` field 
in this [spreadsheet](https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345) 
This sample file is already filtered to those samples having the introgressed teosinte haplotype for this region of the genome.
We made this table using vlooukps in google sheets.

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

We will use this file to filter  `hpc1_aligned.fasta`
```
# results of
# grep '>' hpc1_aligned.fasta
>hpc1_B73  Zm00001eb121780 3:8495640-8499737 dna:chromosome chromosome:Zm-B73-REFERENCE-NAM-5.0:3:1:238017767:1 REF
>hpc1_S29
>hpc1_S82
>hpc1_S26
>hpc1_S41
>hpc1_S72
>hpc1_S45
```
Using the following command:

```bash
# Convert fasta files to one line sequence files
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf "%s",$0}}}' hpc1_aligned.fasta > hpc1_oneline.fasta
```

Additionally we need a table pairing the sequence name in the `hpc1_oneline.fasta` with the
the `label` field from the [spreadsheet](https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345) 
We also made this table using vlooukps in google sheets.

#### `name_swap.tab`
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

Now we use the info in `name_swap.tab` to have visualize the `taxa` field in the alignments and phylogenetic trees

```bash
# Convert fasta files to one line sequence files
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf "%s",$0}}}' hpc1_aligned.fasta > hpc1_oneline.fasta

# Create empty file
echo -n > hpc1_filtered.fasta

# Add reference sequences
grep -A1 -p "B73" hpc1_oneline.fasta  >> hpc1_filtered.fasta
grep -A1 -p "TIL18" hpc1_oneline.fasta  >> hpc1_filtered.fasta

# Filter samples
cut -f1 sample_label.tab | while IFS=$'\t'  read -r SAMPLE;  do grep -A1 -p "${SAMPLE}$" hpc1_oneline.fasta >> hpc1_filtered.fasta; done

# redo alignment
mafft --reorder --adjustdirection --nuc --auto hpc1_filtered.fasta > hpc1_filtered_realigned.fasta

# Swap names in the fasta header
awk 'FNR==NR{a[$1]=$2;next} /^>/{$0=">"a[substr($0,2)];} 1'  name_swap.tab hpc1_filtered_realigned.fasta

# B73 and TIL18 had to be adjusted manually   
```

at the end you should have the folloeing files:

```
name_swap
├── hpc1_aligned.fasta
├── hpc1_filtered.fasta
├── hpc1_filtered_realigned.fasta
├── hpc1_nice_labels.fasta
├── hpc1_oneline.fasta
├── name_swap.tab
└── sample_label.tab
```

### Name swapping command explained:
1. `FNR==NR{a[$1]=$2;next}` - Reads the name_swap.tab file and stores mappings in array `a`
2. `/^>/{$0=">"a[substr($0,2)];}` - For header lines, replaces entire line with new header
3. `substr($0,2)` - Extracts the sequence ID without the ">" character
4. `1` - Prints every line (modified or not)

### Alternative: R Biostrings Approach
For a more readable solution in R, consider using the Biostrings package for FASTA manipulations.
