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
ipt6	Zm00001eb062030
hpc1	Zm00001eb121780
nrg11	Zm00001eb206940
nlp15	Zm00001eb231720
gdsl	Zm00001eb268440
tcptf9	Zm00001eb372490
```

### 2. Create Taxa Reference Table


Create a tab-separated file named `taxa_db.tab` with columns for:
1. Taxa short name (used for file naming)
2. Genome file prefix (used for BLAST database)
3. Gene database name (used for transcript extraction)
4. Column position in orthogroup table (for mapping genes)

```
B73	Zm-B73-REFERENCE-NAM-5.0	Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gene	2
Gigi	Zd-Gigi-REFERENCE-PanAnd-1.0	Zd-Gigi-REFERENCE-PanAnd-1.0_Zd00001aa.1.gene	7
Momo	Zd-Momo-REFERENCE-PanAnd-1.0	Zd-Momo-REFERENCE-PanAnd-1.0_Zd00003aa.1.gene	8
Zh	Zh-RIMHU001-REFERENCE-PanAnd-1.0	Zh-RIMHU001-REFERENCE-PanAnd-1.0_Zh00001aa.1.gene	9
Zn	Zn-PI615697-REFERENCE-PanAnd-1.0	Zn-PI615697-REFERENCE-PanAnd-1.0_Zn00001aa.1.gene	10
TIL01	Zv-TIL01-REFERENCE-PanAnd-1.0	Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.gene	4
TIL11	Zv-TIL11-REFERENCE-PanAnd-1.0	Zv-TIL11-REFERENCE-PanAnd-1.0_Zv00002aa.1.gene	6
TIL18	Zx-TIL18-REFERENCE-PanAnd-1.0	Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002aa.1.gene	3
TIL25	Zx-TIL25-REFERENCE-PanAnd-1.0	Zx-TIL25-REFERENCE-PanAnd-1.0_Zx00003aa.1.gene	5
TdFL	Td-FL_9056069_6-REFERENCE-PanAnd-2.0a	Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_Td00001bc.1.gene	11
CML457	Zm-CML457-REFERENCE-HiLo-1.0	Zm-CML457-REFERENCE-HiLo-1.0_Zm00106aa.1.gene	12
CML459	Zm-CML459-REFERENCE-HiLo-1.0	Zm-CML459-REFERENCE-HiLo-1.0_Zm00107aa.1.gene	13
CML530	Zm-CML530-REFERENCE-HiLo-1.0	Zm-CML530-REFERENCE-HiLo-1.0_Zm00108aa.1.gene	14
PDJ	Zm-PDJ-REFERENCE-HiLo-1.0	Zm-PDJ-REFERENCE-HiLo-1.0_Zm00112aa.1.gene	15
PT	Zm-PT-REFERENCE-HiLo-1.0	Zm-PT-REFERENCE-HiLo-1.0_Zm00109aa.1.gene	16
TAB	Zm-TAB-REFERENCE-HiLo-1.0	Zm-TAB-REFERENCE-HiLo-1.0_Zm00111aa.1.gene	17
ZAP	Zm-ZAP-REFERENCE-HiLo-1.0	Zm-ZAP-REFERENCE-HiLo-1.0_Zm00110aa.1.gene	18
```
At some point we want to have all these taxa:

#### Zd Zea diploperennis      ✓
- Zd-Gigi-REFERENCE-PanAnd-1.0 ✓
- Zd-Momo-REFERENCE-PanAnd-1.0 ✓

#### Zh Zea huehuetenangensis
- Zh-RIMHU001-REFERENCE-PanAnd-1.0 ✓

#### Zl Zea luxurians ✘
- Zl-RIL003-REFERENCE-PanAnd-1.0 ✘ Pending annotation from liftovertools

#### Zn Zea nicaraguensis.         ✓
- Zn-PI615697-REFERENCE-PanAnd-1.0 ✓

#### Zv Zea mays spp. parviglumis  ✓
- Zv-TIL01-REFERENCE-PanAnd-1.0    ✓
- Zv-TIL11-REFERENCE-PanAnd-1.0    ✓

#### Zx Zea mays spp. mexicana.    ✓
- Zx-TIL18-REFERENCE-PanAnd-1.0    ✓
- Zx-TIL25-REFERENCE-PanAnd-1.0    ✓

#### Td Tripsacum dactyloides             
- Td-FL_9056069_6-REFERENCE-PanAnd-2.0 ✓
- Td-KS_B6_1-REFERENCE-PanAnd-2.0      ✘ no blast hits for many genes

#### Zea mays from HiLo project ✓

- Zm-CML457-REFERENCE-HiLo-1.0 ✓
- Zm-CML459-REFERENCE-HiLo-1.0 ✓
- Zm-CML530-REFERENCE-HiLo-1.0 ✓
- Zm-PDJ-REFERENCE-HiLo-1.0    ✓
- Zm-PT-REFERENCE-HiLo-1.0     ✓
- Zm-TAB-REFERENCE-HiLo-1.0    ✓
- Zm-ZAP-REFERENCE-HiLo-1.0    ✓


### 2. Download Reference Genomes and Annotations

***The wget commands work only from the login node***

Create a script called `download_references.sh`:

```bash
#!/bin/bash
# download_references.sh
# Script to download genome and annotation files for multiple maize taxa

# Function to download and decompress a file
download_file() {
  local url=$1
  local output_file=$2
  local file_type=$3
  local final_name=$4
  
  if [ ! -f "$final_name" ]; then
    echo "  Downloading $file_type..."
    wget "$url"
    gunzip "$output_file"
    
    # Rename if needed
    if [ "$output_file" != "$final_name" ] && [ -f "$output_file" ]; then
      mv "$output_file" "$final_name"
    fi
  else
    echo "  $file_type already exists"
  fi
}

# Create ref directory if it doesn't exist
mkdir -p ref
cd ref

# Read the taxa_db.tab file to get the list of genomes to download
while IFS=$'\t' read -r GENOME GENOME_PREFIX GENE_DB COL; do
  echo "Checking files for $GENOME..."
  
  # Extract GFF3 prefix from gene DB name
  GFF3_PREFIX=$(echo $GENE_DB | sed 's/\.gene$//')
  
  # Determine base URL and file paths based on genome
  if [[ "$GENOME" == "TdFL" ]]; then
    BASE_URL="https://download.maizegdb.org/Td-FL_9056069_6-REFERENCE-PanAnd-2.0/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a"
    GENOME_FASTA="${GENOME_PREFIX}.fa.gz"
    GENOME_OUTPUT="${GENOME_PREFIX}.fa"
  elif [[ "$GENOME" == "TdKS" ]]; then
    BASE_URL="https://download.maizegdb.org/Td-KS_B6_1-REFERENCE-PanAnd-2.0/Td-KS_B6_1-REFERENCE-PanAnd-2.0a"
    GENOME_FASTA="Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fa.gz"
    GENOME_OUTPUT="Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fa"
  else
    BASE_URL="https://download.maizegdb.org/${GENOME_PREFIX}"
    GENOME_FASTA="${GENOME_PREFIX}.fa.gz"
    GENOME_OUTPUT="${GENOME_PREFIX}.fa"
  fi
  
  # Download genome FASTA
  download_file "${BASE_URL}/${GENOME_FASTA}" "${GENOME_OUTPUT}" "genome FASTA for $GENOME" "${GENOME_PREFIX}.fa"
  
  # Download annotation GFF3
  download_file "${BASE_URL}/${GFF3_PREFIX}.gff3.gz" "${GFF3_PREFIX}.gff3" "annotation GFF3 for $GENOME" "${GFF3_PREFIX}.gff3"
  
  # Download gene FASTA for ALL genomes
  download_file "${BASE_URL}/${GENE_DB}.fa.gz" "${GENE_DB}.fa" "gene FASTA for $GENOME" "${GENE_DB}.fa"
  
  # PT chromosome numbers are padded in the gff3 but not in the reference
  perl -i -p -e 's/chr0/chr/' Zm-PT-REFERENCE-HiLo-1.0_Zm00109aa.1.gff3
  
done < ../taxa_db.tab

echo "All required files are now available in the ref directory"
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

# Create BLAST databases for all downloaded FASTA files
# These databases enable efficient sequence retrieval

cd ref

# Index all FASTA files for BLAST search and sequence retrieval
for s in *.fa; do
  OUTDB=$(echo $s | sed 's/\.fa$//')
  
  # Check if BLAST index files already exist
  if [ -f "${OUTDB}.nin" ] && [ -f "${OUTDB}.nhr" ] && [ -f "${OUTDB}.nsq" ]; then
    echo "BLAST database for $s already exists, skipping..."
  else
    echo "Creating BLAST database for $s this might take few minutes..."
    makeblastdb -in $s -out $OUTDB -dbtype nucl -parse_seqids
  fi
done

cd ..
```

Make the script executable and run it:

```bash
chmod +x create_blast_dbs.sh
./create_blast_dbs.sh
```

### Why Use BLAST+ Databases?

BLAST+ databases provide several advantages for comparative genomics work:

1. **Efficient sequence storage and indexing**: BLAST databases convert FASTA files into specialized, indexed formats that enable rapid sequence searching and retrieval.

2. **Sequence lookup by ID**: When working with multiple genomes and thousands of genes, direct sequence lookup by identifier is much faster than parsing FASTA files.

3. **Range-based extraction**: BLAST databases allow you to extract specific regions of sequences (like genes with flanking regions) without loading entire chromosomes into memory.

4. **Sequence orientation control**: You can easily extract sequences in either forward or reverse complement orientation based on the strand information.

5. **Pre-formatted for homology searches**: The same databases used for sequence extraction can be used for BLAST searches to identify orthologous sequences.

### How BLAST+ Databases Are Created and Indexed

The `makeblastdb` command is used to create BLAST databases from FASTA files. Here's what happens during the process:

```bash
makeblastdb -in genome.fa -out genome_db -dbtype nucl -parse_seqids
```

- **Database creation**: The command processes the input FASTA file and creates multiple specialized files with extensions like .nhr, .nin, .nsq, which collectively form the BLAST database.

- **Sequence parsing**: The `-parse_seqids` flag is crucial as it instructs makeblastdb to interpret and index the sequence identifiers from the FASTA headers. This enables later retrieval of sequences by their IDs.

- **Indexing mechanism**: The tool creates multiple index files:
  - `.nhr` files contain header information
  - `.nin` files contain the index of sequence names
  - `.nsq` files contain the actual sequence data in a compressed format

- **Database types**: The `-dbtype nucl` parameter specifies that we're working with nucleotide sequences (as opposed to protein sequences with `-dbtype prot`).

Without proper indexing, sequence retrieval would require linear scanning of potentially huge genome files, making the process extremely slow.

### 4. Identify Orthologous Genes
***Work in an interactive session***
We'll need  `taxa_db.tab` again here.

Create a script called `identify_orthologs.sh`:

```bash
#!/bin/bash
# identify_orthologs.sh

# This script identifies orthologous genes across all taxa
# using either OrthoMCL results or BLAST searches

# Make sure B73_gene_targets.tab is tab separated
perl -i -pe 's/ +/\t/' B73_gene_targets.tab

# Extract target B73 gene IDs
cut -f2 B73_gene_targets.tab > B73_gene_targets.list

# Prepare B73 sequences for BLAST
echo "Extracting B73 target sequences..."
blastdbcmd -db ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gene \
  -entry_batch B73_gene_targets.list > B73_gene_targets_genomic.fasta

# Add gene symbols to FASTA headers
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{ sub(/>/,">"a[$1]" ",$1) }1' \
  B73_gene_targets.tab B73_gene_targets_genomic.fasta > B73_gene_targets_genomic.fasta.tmp
mv B73_gene_targets_genomic.fasta.tmp B73_gene_targets_genomic.fasta

# Get gene symbols for the first column (column 01)
cut -f1 B73_gene_targets.tab > c01.tmp

# Use OrthoMCL data for pre-computed orthologs if available
ORTHOMCL_FILE="/rsstu/users/r/rrellan/DOE_CAREER/inv4m/synteny/results/Orthogroups.tsv"
if [ -f "$ORTHOMCL_FILE" ]; then
  echo "Using OrthoMCL results from $ORTHOMCL_FILE"
  
  # Process OrthoMCL data for B73, TIL18, and TIL01 (columns 1-3)
  for t in $(cut -f2 B73_gene_targets.tab); do
    grep $t $ORTHOMCL_FILE | \
      perl -pe 's/,.*?\t/\t/g; s/_P\d+//g' | \
      perl -pe 's/,\s\S+//g' | \
      perl -pe 's/\r\n/\n/g' | \
      cut -f2-4 >> c02.tmp
  done
fi

# Loop through taxa_db.tab and run BLAST against each taxon
# Skip B73, TIL18, and TIL01
while IFS=$'\t' read -r GENOME GENOME_PREFIX GENE_DB COL; do
  # Ensure column number is padded to 2 digits
  COL_PADDED=$(printf "%02d" $COL)
  
  # Skip B73, TIL18, and TIL01
  if [[ "$GENOME" != "B73" && "$GENOME" != "TIL18" && "$GENOME" != "TIL01" ]]; then
    echo "BLASTing against ${GENOME}..."
    
    # Check if the gene database exists
    if [ -f "ref/${GENE_DB}.fa" ]; then
      TARGET_DB="ref/${GENE_DB}"
    else
      # Use genome database if gene database doesn't exist
      TARGET_DB="ref/${GENOME_PREFIX}"
      echo "Note: Using genome database for $GENOME as gene database not found"
    fi
    
    # Run BLAST to find orthologs
    blastn -task megablast \
      -query B73_gene_targets_genomic.fasta \
      -db "$TARGET_DB" \
      -outfmt 6 \
      -num_alignments 1 -max_hsps 1 \
      -num_threads 4 > ${GENOME}.blast

    # Extract best hits (filtering out problematic hits)
    grep -v "gpat" ${GENOME}.blast | \
      grep -v "aaap69" | \
      cut -f2 | \
      sed 's/::.*//' > c${COL_PADDED}.tmp
  else
    echo "Skipping BLAST for ${GENOME} (using OrthoMCL data or already processed)"
  fi
done < taxa_db.tab

# Combine all orthologs into final table
paste c*.tmp > gene_targets_orthogroup.tab

# Clean up temporary files
rm c*.tmp *.blast

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

#### Understanding GFF3 Format

The GFF3 (Generic Feature Format version 3) format is a tab-delimited text file used to describe genomic features. Each feature is described on a single line with 9 fields:

1. **Sequence ID**: Chromosome or scaffold identifier
2. **Source**: Program or database that generated the feature
3. **Type**: Feature type (e.g., gene, mRNA, exon, CDS)
4. **Start**: 1-based start coordinate of the feature
5. **End**: End coordinate of the feature (inclusive)
6. **Score**: Numeric score (or "." if not applicable)
7. **Strand**: "+" for forward strand, "-" for reverse strand
8. **Phase**: For CDS features, indicates where the next codon starts (0, 1, 2, or ".")
9. **Attributes**: Semicolon-separated list of tag=value pairs with additional information

Example GFF3 entry for an mRNA feature:
```
Chr1  maker  mRNA  12345  23456  .  +  .  ID=Zm00001eb062030;Parent=gene1;canonical_transcript=1;Name=ipt6
```

In this pipeline, we're specifically looking for lines with:
- Feature type "mRNA" (field 3)
- The gene ID we're targeting (in the Attributes field)
- The attribute "canonical_transcript=1", which indicates this is the primary transcript for the gene

The canonical transcript is important because many genes have alternative splicing that produces multiple transcripts, but we only want to extract one representative transcript per gene.

***Work in an interactive session***


Create a script called `extract_canonical_transcripts.sh`:

```bash
#!/bin/bash
# extract_canonical_transcripts.sh

# This script extracts canonical transcript coordinates for each gene
# across all taxa from their respective GFF3 annotation files

# Read taxa information from taxa_db.tab
while IFS=$'\t' read -r GENOME GENOME_PREFIX GENE_DB COL; do
  echo "Extracting ${GENOME} canonical transcripts..."

  # Extract GFF3 prefix from gene DB name (remove .gene suffix)
  GFF3_PREFIX=$(echo $GENE_DB | sed 's/\.gene$//')
  echo "${GFF3_PREFIX}.gff3"

  # Create empty GFF3 file for this taxon
  > ${GENOME}_gene_targets.gff3

  # Extract canonical transcripts for each gene
  # Use the appropriate column from orthogroup table based on the COL value
  # Ensure column number is padded to 2 digits for consistency
  COL_PADDED=$(printf "%02d" $COL)

  for t in $(cut -f${COL} gene_targets_orthogroup.tab); do
    # Skip empty entries
    if [[ "$t" != "." ]]; then
      # Find canonical transcript and append to GFF3 file
      # Use the GFF3_PREFIX derived from GENE_DB column for the gff3 file
      grep -w $t ref/${GFF3_PREFIX}.gff3 | \
        grep mRNA | \
        grep "canonical_transcript=1" >> ${GENOME}_gene_targets.gff3
    fi
  done

  # Report number of transcripts found
  COUNT=$(wc -l < ${GENOME}_gene_targets.gff3)
  echo "  Found $COUNT canonical transcripts for $GENOME"

done < taxa_db.tab

echo "Canonical transcript annotations in GFF3 files"

# Show summary of GFF files
for GENOME in $(cut -f1 taxa_db.tab); do
  COUNT=$(wc -l < ${GENOME}_gene_targets.gff3)
  echo "$GENOME: $COUNT transcripts"
done
```

### 6. Create Coordinate Files with Flanking Regions

#### How blastdbcmd Retrieves Subsequences

The `blastdbcmd` tool is a powerful utility for extracting sequences or subsequences from BLAST databases. 
It used in the following manner:
```bash
blastdbcmd -db ref/Zm-B73-REFERENCE-NAM-5.0 -entry_batch B73_gene_targets_entry_batch.tab > B73_targets.fas
```

We will use the gff files from the previous step to build `entry_batch` files for sequence retrieval.
Each `entry_batch` file has the necessary info per line to retrieve a specific subsequence from a blast database.

Here's how `blastdbcmd` works:

1. **Sequence retrieval mechanism**: 
   - blastdbcmd uses the indices created by makeblastdb to directly access sequence data
   - It can retrieve entire sequences or specific regions (subsequences) based on coordinates

2. **The entry_batch format**:
   - Each line contains a sequence identifier followed by optional range and strand specifications
   - Format: `seqid range strand`
   - Example: `chromosome01 1000-2000 plus`
   - This means "extract bases 1000 through 2000 from chromosome01 in the forward orientation"

   https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastdbcmd_application_opti/

3. **Coordinate system**:
   - BLAST databases use 1-based inclusive coordinates (first base is position 1)
   - The range specification is formatted as `start-end` 
   - For genes on the minus strand, you specify `minus` to get the reverse complement

4. **Handling subsequence extraction**:
   - When we extract genes with 2kb flanking regions, blastdbcmd performs several steps:
     - Locates the sequence (e.g., chromosome) containing the gene
     - Extracts the specified range (gene ± 2kb)
     - If on the minus strand, automatically reverse-complements the sequence
     - Returns the subsequence with a header containing the original sequence ID and coordinates

5. **Efficiency benefits**:
   - Direct byte-offset addressing allows retrieval of just the relevant portion of a chromosome
   - No need to load entire genome sequences into memory
   - Supports batch processing of multiple sequence requests at once

The entry_batch approach is particularly valuable in this pipeline because we're extracting multiple gene sequences from multiple reference genomes, each with specific coordinates and orientations. Without blastdbcmd, we would need to write custom parsers to handle sequence extraction from FASTA files, which would be more error-prone and less efficient.


***Work in an interactive session***

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
See section 6 for full explanation of `blastdbcmd` fro sequence retrieval.

***Work in an interactive session***


Then create the main script to extract sequences:

```bash
#!/bin/bash
# get_target_sequences.sh
# this is a slow step maybe do it with samtools

REF_DIR="./ref"
TAXA_DB="taxa_db.tab"

# Process each taxon
while IFS=$'\t' read -r TAXA DB GENEDB COL; do
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

bash create_blast_dbs.sh

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

3. **Create BLAST databases**: We index the genomes to enable BLAST searches and efficient sequence retrieval. The BLAST+ database format provides random access to sequences without needing to parse entire FASTA files.

4. **Identify orthologous genes**: We use either:
   - Pre-computed OrthoMCL results if available
   - Direct BLAST searches to find the best hit in each species

5. **Extract canonical transcripts**: For each gene in each species, we extract the coordinates of the canonical transcript from the GFF3 annotation file.

6. **Add flanking regions**: We extend the coordinates by 2kb upstream and 2kb downstream to capture potential regulatory regions.

7. **Extract target sequences**: We use blastdbcmd to extract the sequences based on the extended coordinates. This leverages the indexed database structure to efficiently retrieve specific subsequences from potentially very large genome files.

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
  - Verify the BLAST database files (.nhr, .nin, .nsq) exist and are not corrupted

- **Empty GFF files**: 
  - Verify that orthologs exist in all species
  - Some genes may not have clear orthologs in more distant species

- **Permission denied errors**: 
  - Make scripts executable with `chmod +x *.sh`

- **Windows line endings**:
  - If scripts fail with syntax errors, check for Windows line endings
  - Fix with `dos2unix *.sh`

- **blastdbcmd errors**:
  - If you get "Entry not found" errors, verify the sequence IDs exist in the original FASTA
  - Make sure reference databases were created with the `-parse_seqids` flag
  - Check that chromosome/scaffold IDs match exactly between your entry_batch file and the database

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

- **BLAST database optimization**:
  - For very large genomes, use the `-max_file_sz` parameter with makeblastdb to control database file size
  - Use `-hash_index` for larger genomes when memory is limited during database creation

## Further Analysis Ideas

- Multiple sequence alignment of orthologous genes
- Phylogenetic analysis of gene families
- Identification of conserved non-coding sequences
- Promoter analysis for transcription factor binding sites
- Structural variation detection between species
