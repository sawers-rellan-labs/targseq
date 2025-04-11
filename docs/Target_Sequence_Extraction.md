# Maize Comparative Genomics: Target Sequence Extraction Pipeline

## Overview

This tutorial guides you through extracting homologous gene sequences (including flanking regions) from multiple Zea taxa and a related Tripsacum species. Starting with a set of B73 maize gene identifiers, you'll:

1. Download reference genomes and annotations
2. Create BLAST databases
3. Identify orthologous genes across taxa
4. Extract target sequences with 2kb flanking regions
5. Rename sequences for easier interpretation

### 0 Setup

#### Login to hazel
```
ssh user@login.hpc.ncsu.edu
```
Downloading and saving files must be done in the login node.

***Please use an interactive session after step 2***

## Prerequisites

- Basic command line skills
- Installed programs:
  - BLAST+ suite (makeblastdb, blastn, blastdbcmd)
  - Standard Unix tools (awk, grep, cut, paste, sed, perl)

## Workflow Steps

### 1. Prepare Your Gene Target List

Create a tab-separated file named `B73_gene_targets.tab` with B73 gene IDs and gene symbols:

```
Zm00001eb062030 ipt6
Zm00001eb121780 hpc1
Zm00001eb206940 nrg11
Zm00001eb231720 nlp1
Zm00001eb268440 gdsl
Zm00001eb372490 tcptf9
```

### 2. Download Reference Genomes and Annotations

***The wget commands work only from the login node***

Create a script called `download_references.sh`:

```bash
#!/bin/bash
# download_references.sh

# Create working directory
mkdir -p ref
cd ref

# B73 maize reference genome
wget 'https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz'
wget 'https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gene.fa.gz'
wget 'https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz'

# Teosinte parviglumis (TIL01)
wget 'https://download.maizegdb.org/Zv-TIL01-REFERENCE-PanAnd-1.0/Zv-TIL01-REFERENCE-PanAnd-1.0.fa.gz'
wget 'https://download.maizegdb.org/Zv-TIL01-REFERENCE-PanAnd-1.0/Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.gff3.gz'

# Teosinte mexicana (TIL18)
wget 'https://download.maizegdb.org/Zx-TIL18-REFERENCE-PanAnd-1.0/Zx-TIL18-REFERENCE-PanAnd-1.0.fa.gz'
wget 'https://download.maizegdb.org/Zx-TIL18-REFERENCE-PanAnd-1.0/Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.gff3.gz'

# Tripsacum dactyloides (Florida accession, because Z)
wget 'https://download.maizegdb.org/Td-FL_9056069_6-REFERENCE-PanAnd-2.0/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fa.gz'
wget 'https://download.maizegdb.org/Td-FL_9056069_6-REFERENCE-PanAnd-2.0/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_Td00001bc.1.gff3.gz'
wget 'https://download.maizegdb.org/Td-FL_9056069_6-REFERENCE-PanAnd-2.0/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_Td00001bc.1.gene.fa.gz'

# Uncompress all files
echo "Uncompressing files, this might take a couple of minutes..."
gunzip *.gz
echo "Done!"

cd ..
```

Make the script executable and run it:

```bash
chmod +x download_references.sh
./download_references.sh
```




### 3. Create BLAST Databases
***From this point onwards use an interactive session***

#### Request an interactive session
Although this pipeline doesn't take much resources it's better for you to avoid working in the login node and use dedicated cluster machines.

```
bsub -Is -n 8 -R "rusage[mem=16]" -W 8:00 bash
```

#### Activate `secapr` conda environment
You will be using BLAST+ which is installed in the secapr environment

```
conda activate /share/maize/user/conda/env/secapr_env
```

Create a script called `create_blast_dbs.sh`:

```bash
#!/bin/bash
# create_blast_dbs.sh

cd ref

# Index all FASTA files for BLAST search and sequence retrieval
for s in *.fa; do
  OUTDB=$(echo $s | sed 's/\.fa$//')
  echo "Creating BLAST database for $s..."
  makeblastdb -in $s -out $OUTDB -dbtype nucl -parse_seqids
done

cd ..
```

Make the script executable and run it:

```bash
chmod +x create_blast_dbs.sh
./create_blast_dbs.sh
```

### 4. Identify Orthologous Genes
***Work in an interactive session***

Create a script called `identify_orthologs.sh`:

```bash
#!/bin/bash
# identify_orthologs.sh

# Mahe sure  B73_gene_targets.tab is tab separated
perl -i -pe 's/ +/\t/' B73_gene_targets.tab

# Extract target B73 gene IDs
cut -f1 B73_gene_targets.tab > B73_gene_targets.list

# Prepare B73 sequences for BLAST
echo "Extracting B73 target sequences..."
blastdbcmd -db ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gene \
  -entry_batch B73_gene_targets.list > B73_gene_targets_genomic.fasta.tmp

# Add gene symbols to FASTA headers
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{ sub(/>/,">"a[$1]" ",$1) }1' \
  B73_gene_targets.tab B73_gene_targets_genomic.fasta.tmp > B73_gene_targets_genomic.fasta

ORTHOMCL_FILE="/rsstu/users/r/rrellan/DOE_CAREER/inv4m/synteny/results/Orthogroups.tsv"

echo "Using OrthoMCL results from $ORTHOMCL_FILE"
  
  # Create temporary files for collecting ortholog info
  > c2.tmp
  
  # For each target gene, find its orthologs in the OrthoMCL file
  # select just one ortholog, remove protein suffix (_P001)
  # remove problematic spaces
  # convert from windows line endings to unix line endings
  for t in $(cut -f1 B73_gene_targets.tab); do
    grep $t $ORTHOMCL_FILE | \
      perl -pe 's/,.*?\t/\t/g; s/_P\d+//g' | \
      perl -pe 's/,\s\S+//g' | \
      perl -pe 's/\r\n/\n/g' | \
      cut -f2-4 >> c2.tmp
  done
  
  # Get gene symbols for the first column
  cut -f2 B73_gene_targets.tab > c1.tmp

# BLAST against Tripsacum 
# Omit when you run OrthoMCL on all taxa including Tripsacum

echo "BLASTing against Tripsacum..."
blastn -task megablast \
  -query B73_gene_targets_genomic.fasta \
  -db ref/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_Td00001bc.1.gene \
  -outfmt 6 \
  -num_alignments 1 -max_hsps 1 \
  -num_threads 4 > TdFL.blast

# Extract Tripsacum best hits (filtering out problematic hits)
grep -v "gpat" TdFL.blast | grep -v "aaap69" | cut -f2 | sed 's/::.*//' > c3.tmp

# Combine all orthologs into final table
paste c1.tmp c2.tmp c3.tmp > gene_targets_orthogroup.tab

# Clean up temporary files
rm *.tmp

echo "The ortholog table was saved to gene_targets_orthogroup.tab:"
echo ""
head gene_targets_orthogroup.tab
```

Make the script executable and run it:

```bash
chmod +x identify_orthologs.sh
./identify_orthologs.sh
```

### 5. Extract Canonical Transcripts from Annotations
***Work in an interactive session***

Create a script called `extract_canonical_transcripts.sh`:

```bash
#!/bin/bash
# extract_canonical_transcripts.sh

# Extract B73 canonical transcripts
echo "Extracting B73 canonical transcripts..."

# Create empty gff3
> B73_gene_targets.gff3

for t in $(cut -f2 gene_targets_orthogroup.tab); do
  grep -w $t ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | \
  grep mRNA | \
  grep "canonical_transcript=1" >> B73_gene_targets.gff3
done

# Extract mexicana canonical transcripts
echo "Extracting TIL18 canonical transcripts..."
# Create empty gff3
> TIL18_gene_targets.gff3

for t in $(cut -f3 gene_targets_orthogroup.tab); do
  grep -w $t ref/Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.gff3 | \
  grep mRNA | \
  grep "canonical_transcript=1" >> TIL18_gene_targets.gff3
done

# Extract parviglumis canonical transcripts
echo "Extracting TIL01 canonical transcripts..."

# Create empty gff3
> TIL01_gene_targets.gff3

for t in $(cut -f4 gene_targets_orthogroup.tab); do
  grep -w $t ref/Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.gff3 | \
  grep mRNA | \
  grep "canonical_transcript=1" >> TIL01_gene_targets.gff3
done

# Extract Tripsacum canonical transcripts
echo "Extracting TdFL canonical transcripts..."

# Create empty gff3
> TdFL_gene_targets.gff3

for t in $(cut -f5 gene_targets_orthogroup.tab); do
  grep -w $t ref/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_Td00001bc.1.gff3 | \
  grep mRNA | \
  grep "canonical_transcript=1" >> TdFL_gene_targets.gff3
done

echo "Canonical transcript annotations in gff3 files"

# Show file name and first line of the gff files
head -n 1 *_gene_targets.gff3
```

Make the script executable and run it:

```bash
chmod +x extract_canonical_transcripts.sh
./extract_canonical_transcripts.sh
```

### 6. Create Coordinate Files with Flanking Regions
***Work in an interactive session***

We will use the gff files from the previous step to build entry batch files for sequence retrieval with `blastdbcmd`.
Each entry batch file has the necessary info to retrieve a specific subsequence from a blast database.

<blockquote>

entry_batch

Input file for batch processing. The format requires one entry per line; each line should begin with the sequence ID followed by any of the following optional specifiers (in any order): range (format: ‘from-to’, inclusive in 1-offsets), strand (‘plus’ or ‘minus’), or masking algorithm ID (integer value representing the available masking algorithm). Omitting the ending range (e.g.: ‘10-‘) is supported, but there should not be any spaces around the ‘-‘.
</blockquote>

See https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastdbcmd_application_opti/

Create a script called `extract_coordinates.sh`:

```bash
#!/bin/bash
# extract_coordinates.sh

# Add 2kb flanking regions to coordinates for each species
for GENOME in B73 TIL01 TIL18 TdFL; do
  echo "Extracting coordinates with flanking regions for $GENOME..."
  
  # Create coordinates with 2kb upstream and 2kb downstream
  awk '{ 
    $4 = $4 - 2000; 
    if ($4 < 1) $4 = 1;  # Ensure start position is not negative
    $5 = $5 + 2000; 
    print $1,$4 "<_>" $5,$7
  }' ${GENOME}_gene_targets.gff3 | \
  perl -pe 's/\+/plus/; s/-/minus/' | \
  perl -pe 's/\<_>/-/' \
  > ${GENOME}_gene_targets_entry_batch.tab
done

# Show file name and first line of the entry_batch files
head -n 1 *_gene_targets_entry_batch.tab
```

Make the script executable and run it:

```bash
chmod +x extract_coordinates.sh
./extract_coordinates.sh
```

### 7. Extract Target Sequences
***Work in an interactive session***

First, create a tab-separated file named `taxa_db.tab` with columns for the taxa name, the database name, and its correspnding `gene_targets_orthogroup.tab` column index.

```
B73	Zm-B73-REFERENCE-NAM-5.0	2
TIL18	Zx-TIL18-REFERENCE-PanAnd-1.0	3
TIL01	Zv-TIL01-REFERENCE-PanAnd-1.0	4
TdFL	Td-FL_9056069_6-REFERENCE-PanAnd-2.0a	5
```


Then create the main script to extract sequences:

```bash
#!/bin/bash
# get_target_sequences.sh

REF_DIR="./ref"
TAXA_DB="taxa_db.tab"

# Process each taxon
while IFS=$'\t' read -r TAXA DB COL; do
  echo "Retrieving sequences from ${DB}..."
  
  # Extract sequences using coordinates
  blastdbcmd -db ${REF_DIR}/${DB} \
    -entry_batch ${TAXA}_gene_targets_entry_batch.tab > ${TAXA}_targets.fas.tmp
  
  # Get coordinate IDs from sequence headers
  grep '>' ${TAXA}_targets.fas.tmp | perl -pe 's/>| .*//g' > c1.tmp
  
  # Get gene IDs for this taxon from orthogroup file
  cut -f${COL} gene_targets_orthogroup.tab > c2.tmp
  
  # Create mapping from coordinates to gene IDs
  paste c1.tmp c2.tmp > coord_to_geneid.tmp
  
  # Replace sequence headers with gene IDs
  awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{ sub(/>/,">"a[$1]" ",$1) }1' \
    coord_to_geneid.tmp ${TAXA}_targets.fas.tmp > with_geneid.fas.tmp
  
  # Get gene symbols (first column of orthogroup table)
  cut -f1 gene_targets_orthogroup.tab > c2.tmp
  
  # Get gene IDs for this taxon again
  cut -f${COL} gene_targets_orthogroup.tab > c1.tmp
  
  # Create mapping from gene IDs to symbols
  paste c1.tmp c2.tmp > geneid_to_symbol.tmp
  
  # Add gene symbols to sequence headers
  awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{ sub(/>/,">"a[$1]" ",$1) }1' \
    geneid_to_symbol.tmp with_geneid.fas.tmp > ${TAXA}_target_sequences.fasta
  
  # Add taxon identifier to sequence headers
  perl -i -pe "if(/>/){s/ /_${TAXA} /}" ${TAXA}_target_sequences.fasta
  
  # Clean up temporary files
  rm *.tmp
  
done < ${TAXA_DB}

head -n1 *_target_sequences.fasta 

echo "Target sequence extraction complete!"
```

Make the script executable and run it:

```bash
chmod +x get_target_sequences.sh
./get_target_sequences.sh
```

### 8. Complete Pipeline Script

Here's a master script that runs all steps in sequence:

```bash
#!/bin/bash
# run_sequence_extraction_pipeline.sh

echo "1. Creating BLAST databases..."
# bash create_blast_dbs.sh

echo "2. Identifying orthologous genes..."
bash identify_orthologs.sh

echo "3. Extracting canonical transcripts..."
bash extract_canonical_transcripts.sh

echo "4. Creating coordinate files with flanking regions..."
bash extract_coordinates.sh

echo "5. Extracting target sequences..."
bash get_target_sequences.sh

echo "Pipeline completed successfully!"
echo "Output: *_target_sequences.fasta files for each taxon"
ls -1 *_target_sequences.fasta
```

Make the script executable and run it:

```bash
chmod +x run_sequence_extraction_pipeline.sh
./run_sequence_extraction_pipeline.sh
```

## Understanding the Workflow

This pipeline performs the following steps:

1. **Start with a list of B73 maize genes**: We begin with a list of B73 gene IDs and their corresponding gene symbols.

2. **Download reference genomes**: We download genome assemblies and annotations for B73 maize, two teosinte species (parviglumis and mexicana), and Tripsacum.

3. **Create BLAST databases**: We index the genomes to enable BLAST searches.

4. **Identify orthologous genes**: We use either:
   - Pre-computed OrthoMCL results if available
   - Direct BLAST searches to find the best hit in each species

5. **Extract canonical transcripts**: For each gene in each species, we extract the coordinates of the canonical transcript from the GFF3 annotation file.

6. **Add flanking regions**: We extend the coordinates by 2kb upstream and 2kb downstream to capture potential regulatory regions.

7. **Extract target sequences**: We use blastdbcmd to extract the sequences based on the extended coordinates.

8. **Rename sequences for clarity**: We modify the sequence headers to include:
   - The gene symbol from B73
   - The species identifier
   - The gene ID in that species

The final output consists of FASTA files for each species containing the target gene sequences with 2kb flanking regions. These sequences can be used for comparative genomics analyses, such as examining conservation of regulatory elements or evolutionary changes in gene structure.

## Example Output Format

```
>ipt6_B73 Zm00001eb062030 
ACTGGCATTCGCTAGCTAGCTGA...
>ipt6_TIL01 Zv00001aa123456
ACTGGCATTCGCTAGCTAGCTGA...
```

## Troubleshooting

- **Missing sequences**: 
  - Check that gene IDs match between your target list and annotation files
  - Verify that the genes have canonical transcripts annotated

- **BLAST database errors**: 
  - Ensure BLAST databases have been created with `-parse_seqids` flag
  - Check that the paths to reference files are correct

- **Empty GFF files**: 
  - Verify that orthologs exist in all species
  - Some genes may not have clear orthologs in more distant species

- **Permission denied errors**: 
  - Make scripts executable with `chmod +x *.sh`

- **Windows line endings**:
  - If scripts fail with syntax errors, check for Windows line endings
  - Fix with `dos2unix *.sh`

## Advanced Options

- **Modify flanking region size**:
  - Edit the `extract_coordinates.sh` script to change the 2kb value
  - For promoter-only analysis, use asymmetric flanking (e.g., 2kb upstream, 0kb downstream)

- **Add more species**:
  - Update the `download_references.sh` script to download additional genomes
  - Add the new species to the `taxa_db.tab` file
  - Repeat the ortholog identification step for the new species

- **Custom filtering**:
  - Add filters in the `identify_orthologs.sh` script to select hits based on criteria like identity percentage
  - Modify the `extract_canonical_transcripts.sh` script to select specific transcript variants

- **Parallel processing**:
  - Add GNU Parallel support for faster processing of multiple genes
  - Example: `cat gene_list | parallel -j 8 './process_gene.sh {}'`

## Further Analysis Ideas

- Multiple sequence alignment of orthologous genes
- Phylogenetic analysis of gene families
- Identification of conserved non-coding sequences
- Promoter analysis for transcription factor binding sites
- Structural variation detection between species