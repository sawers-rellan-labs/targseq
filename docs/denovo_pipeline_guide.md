# Step by step guide for *de novo* assembly of target sequencing data

This guide explains a bioinformatics pipeline designed for processing Illumina sequencing data, performing de novo assembly, identifying relevant sequences through BLAST searches, and creating multiple sequence alignments. The pipeline uses LSF (Load Sharing Facility) for job management in a high-performance computing environment.

## Set up
We installed `secapr` from:
https://github.com/AntonelliLab/seqcap_processor

See `secapr_installation.md`

#### Login to hazel
```
ssh user@login.hpc.ncsu.edu
```
As we are sending jobs in batch through LSF `bsub` to be run in the cluster,
you should execute this pipeline from the login node.

#### Execute code from your user folder in targseq
```
cd /rsstu/users/r/rrellan/BZea/targseq/
cd user
```
#### Copy example data to your user folder in targseq

```
mkdir raw
ls -1 ../data | grep -P "_S20_|_S27_|_S29_" | xargs -I{} cp ../data/{} raw
```

## Overview of the Pipeline

The pipeline consists of four main scripts that run sequentially:

1. **q_clean_reads.sh**: Cleans raw Illumina sequencing reads using fastp
2. **q_assemble_reads.sh**: Performs de novo assembly of the cleaned reads using SPAdes
3. **q_find_target_contigs.sh**: Conducts BLAST searches to identify sequences of interest
4. **q_align_sequences.sh**: Creates multiple sequence alignments using MAFFT

Each script builds upon the output of the previous script, and they should be executed in order.
The prefix `q_` stands for queuing LSF jobs. 
These batch scripts are wrappers for the job submission command `bsub`.

```bash
# Add executable permission to shell scripts
chmod +x batch_scripts/q_*.sh
./batch_scripts/q_clean_reads.sh
# wait for the jobs to run then do:
./batch_scripts/q_assemble_reads.sh
# wait for the jobs to run then do:
./batch_scripts/q_find_target_contigs.sh
# wait for the jobs to run then do:
./batch_scripts/q_align_sequences.sh
```

## Prerequisites

Before starting, ensure you have:

- A conda environment with the necessary tools installed (secapr_env in this example)
- Raw Illumina sequencing data (paired-end reads in fastq.gz format)
- A reference sequences file (B73_target_sequences.fasta)
- A sample annotation file (sample_annotation.tab) containing sample IDs and filenames
- Access to an LSF job scheduling system

### Reference sequence files

The file  `B73_target_sequences.fasta` comes from the B73 genome annotation `Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3`
found in https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/.

See `Target_Sequence_Extraction.md`.

The sequences include the canonical transcript sequence plus 2Kb upstream the TSS as described by the following coordinates:

```
# canonical transcript gene body coordinates extracted from Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
# chr start end strand gene symbol
1 298566258 298567661 plus Zm00001eb062030 ipt6
3 8497640 8499737 plus Zm00001eb121780 hpc1
4 244043231 244045806 plus Zm00001eb206940 nrg11
5 78778520 78781733 minus Zm00001eb231720 nlp1
6 64573577 64575876 plus Zm00001eb268440 gdsl
9 6311753 6313844 minus Zm00001eb372490 tcptf9

# 2kb 'promoter' + gene body coordinates for blastdbcmd
1 298564258-298567661 plus
3 8495640-8499737 plus
4 244041231-244045806 plus
5 78778520-78783733 minus
6 64571577-64575876 plus
9 6311753-6315844 minus
```

Similarly the sequences from TIL18, an inbred line of *Zea mays spp. mexicana*, 
are stored in  `TIL18_target_sequences.fasta` and were extracted from  `Zx-TIL18-REFERENCE-PanAnd-1.0_Zx00002ab.1.gff3`
in https://download.maizegdb.org/Zx-TIL18-REFERENCE-PanAnd-1.0/

See `Target_Sequence_Extraction.md`

### Sample annotation file

The sample annotation file should be tab-separated with at least two columns: Sample ID and filename. Example:

```
S20	Zx0550_P4_P3_P5311_S20_L001_R1_001.fastq.gz
S20	Zx0550_P4_P3_P5311_S20_L001_R2_001.fastq.gz
S27	Zx0540_P3_P5_P1111_S27_L001_R1_001.fastq.gz
S27	Zx0540_P3_P5_P1111_S27_L001_R2_001.fastq.gz
S29	Zx0580_P2_P5_P2411_S29_L001_R1_001.fastq.gz
S29	Zx0580_P2_P5_P2411_S29_L001_R2_001.fastq.gz
```

You can make it as a spreadsheet and save it as tsv. 
But we'll run a shell command using perl and regex to generate the annotation file from the folder contents.

```
ls -1 raw/ | perl -pe 's/(^.*(S\d+))_L001/$2\t$1\_L001/' | sort -n -k1.2 > sample_annotation.tab
```

## Directory Structure

The pipeline uses the following directory structure:

```
├── raw/                   # Raw sequencing data
├── clean/                 # Cleaned reads
├── assemblies/            # De novo assemblies
├── blast_results/         # BLAST search results
├── alignments/            # Multiple sequence alignments
├── logs/                  # Log files
├── B73_target_sequences.fasta  # B73 target sequences
├── sample_annotation.tab      # Sample annotations
└── batch_scripts
    ├── q_clean_reads.sh
    ├── q_assemble_reads.sh
    ├── q_find_target_contigs.sh
    └── q_align_sequences.sh
```

## Step 1: Read Cleaning (q_clean_reads.sh)

This script cleans raw Illumina reads using fastp, which performs quality filtering, adapter trimming, and read deduplication.

### Script Explanation

```bash
#!/bin/bash
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/secapr_env

# Define directories and input files
RAW_DIR="raw"
CLEAN_DIR="clean"
ASSEMBLY_DIR="assemblies"
BLAST_DIR="blast_results"
REF_FILE="B73_target_sequences.fasta"
FINAL_HITS="final_combined_hits.fasta"
SAMPLE_FILE="sample_annotation.tab"

# Create output directories if they don't exist
mkdir -p ${CLEAN_DIR}
mkdir -p ${ASSEMBLY_DIR}
mkdir -p ${BLAST_DIR}
mkdir -p logs

# Check if sample file exists
if [ ! -f "${SAMPLE_FILE}" ]; then
    echo "Error: Sample annotation file ${SAMPLE_FILE} not found!"
    exit 1
fi

# Step 1: Clean reads with fastp
echo "Submitting read cleaning jobs..."
while IFS=$'\t' read -r SAMPLE_ID FILENAME; do
    # Extract base filename (without _R1_001.fastq.gz or _R2_001.fastq.gz)
    BASE_FILENAME=$(echo ${FILENAME} | sed 's/_R[12]_001\.fastq\.gz$//')

    # Define input and output files
    R1="${RAW_DIR}/${BASE_FILENAME}_R1_001.fastq.gz"
    R2="${RAW_DIR}/${BASE_FILENAME}_R2_001.fastq.gz"
    OUT_R1="${CLEAN_DIR}/${BASE_FILENAME}_R1_001.clean.fastq.gz"
    OUT_R2="${CLEAN_DIR}/${BASE_FILENAME}_R2_001.clean.fastq.gz"

    # Check if both files exist before submitting job
    if [ -f "${R1}" ] && [ -f "${R2}" ]; then
        # Submit fastp job
        bsub -q "sara" \
             -J "fastp_${SAMPLE_ID}" \
             -o "logs/fastp_${SAMPLE_ID}.out" \
             -e "logs/fastp_${SAMPLE_ID}.err" \
             -n 4 -R "rusage[mem=4G]" \
             -W 60 \
             "fastp -w 4 -i ${R1} -I ${R2} -o ${OUT_R1} -O ${OUT_R2} \
                    --dedup \
                    --qualified_quality_phred 30 \
                    --length_required 50 \
                    --html ${CLEAN_DIR}/${SAMPLE_ID}_fastp.html \
                    --json ${CLEAN_DIR}/${SAMPLE_ID}_fastp.json"

        echo "Submitted fastp job for ${SAMPLE_ID} (${BASE_FILENAME})"
    else
        echo "Warning: Input files for ${SAMPLE_ID} not found. Skipping."
    fi
# Only process lines that have _R1_ in them to avoid duplicates (since file contains both R1 and R2)
done < <(grep "_R1_" ${SAMPLE_FILE})
```

### Running the Script

1. Make sure your raw sequencing data is in the `raw/` directory
2. Ensure your sample annotation file `sample_annotation.tab` is correctly formatted
3. Make the script executable and run it:

```bash
chmod +x ./batch_scripts/q_clean_reads.sh
./batch_scripts/q_clean_reads.sh
```

### What Happens?

- The script reads the sample annotation file and processes each sample
- For each sample, it submits an LSF job that runs fastp with these parameters:
  - Removes duplicate reads (`--dedup`)
  - Requires base quality of 30 or higher (`--qualified_quality_phred 30`)
  - Requires reads to be at least 50 bp long (`--length_required 50`)
  - Uses 4 threads (`-w 4`)
- Output includes cleaned fastq files and HTML/JSON reports
- Jobs are given 60 minutes to complete (`-W 60`)
- Each job uses 4 CPUs and 4GB memory

### Key Parameters

- `-q "sara"`: Specifies the LSF queue to use
- `-J "fastp_${SAMPLE_ID}"`: Names the job with the sample ID for tracking
- `-n 4`: Requests 4 CPUs
- `-R "rusage[mem=4G]"`: Requests 4GB of memory
- `-W 60`: Sets time limit to 60 minutes

## Step 2: De Novo Assembly (q_assemble_reads.sh)

This script performs de novo assembly of the cleaned reads using SPAdes.

### Script Explanation

```bash
#!/bin/bash
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/secapr_env

# Define directories and input files
RAW_DIR="raw"
CLEAN_DIR="clean"
ASSEMBLY_DIR="assemblies"
BLAST_DIR="blast_results"
REF_FILE="B73_target_sequences.fasta"  
FINAL_HITS="final_combined_hits.fasta"
SAMPLE_FILE="sample_annotation.tab"  

# Create output directories if they don't exist
mkdir -p ${CLEAN_DIR}
mkdir -p ${ASSEMBLY_DIR}
mkdir -p ${BLAST_DIR}
mkdir -p logs

# Check if sample file exists
if [ ! -f "${SAMPLE_FILE}" ]; then
    echo "Error: Sample annotation file ${SAMPLE_FILE} not found!"
    exit 1
fi

# Step 2: De novo assembly with SPAdes
echo "Submitting SPAdes assembly jobs..."
while IFS=$'\t' read -r SAMPLE_ID FILENAME; do
    # Extract base filename (without _R1_001.fastq.gz)
    BASE_FILENAME=$(echo ${FILENAME} | sed 's/_R1_001\.fastq\.gz$//')

    # Define input and output files/directories
    CLEAN_R1="${CLEAN_DIR}/${BASE_FILENAME}_R1_001.clean.fastq.gz"
    CLEAN_R2="${CLEAN_DIR}/${BASE_FILENAME}_R2_001.clean.fastq.gz"
    OUT_DIR="${ASSEMBLY_DIR}/${SAMPLE_ID}"

    # Check if cleaned files exist before submitting job
    if [ -f "${CLEAN_R1}" ] && [ -f "${CLEAN_R2}" ]; then
        # Submit SPAdes job with dependency on fastp completion
        bsub -q "sara" \
             -J "spades_${SAMPLE_ID}" \
             -o "logs/spades_${SAMPLE_ID}.out" \
             -e "logs/spades_${SAMPLE_ID}.err" \
             -n 8 -R "rusage[mem=20G]" \
             -W 48:00 \
             "spades.py --isolate -t 8 -m 20 \
                      --only-assembler \
                      --cov-cutoff auto \
                      -1 ${CLEAN_R1} \
                      -2 ${CLEAN_R2} \
                      -o ${OUT_DIR}"

        echo "Submitted SPAdes job for ${SAMPLE_ID}"
    else
        echo "Warning: Cleaned files for ${SAMPLE_ID} not found. Skipping assembly."
    fi
done < <(grep "_R1_" ${SAMPLE_FILE})
```

### Running the Script

Ensure Step 1 has completed, then run:

```bash
chmod +x ./batch_scripts/q_assemble_reads.sh
./batch_scripts/q_assemble_reads.sh
```

### What Happens?

- The script reads the sample file and processes each sample
- It checks for the existence of cleaned read files from Step 1
- For each sample, it submits an LSF job that runs SPAdes with these parameters:
  - Uses the `--isolate` mode (optimized for single isolate assemblies)
  - Sets coverage cutoff to automatic (`--cov-cutoff auto`)
  - Uses 8 threads and 20GB memory
  - Skips read error correction (`--only-assembler`)
- The assembly output is stored in `assemblies/SAMPLE_ID/`
- The main output file of interest is `scaffolds.fasta`

### Key Parameters

- `-n 8 -R "rusage[mem=20G]"`: Requests 8 CPUs and 20GB memory
- `-W 48:00`: Sets a 48-hour time limit
- `--isolate`: Optimizes SPAdes for bacterial/fungal isolate assemblies
- `--only-assembler`: Skips read error correction, which may have already been done in cleaning

## Step 3: BLAST Search (q_find_target_contigs.sh)

This script uses BLAST to identify sequences of interest by comparing the assemblies to reference sequences.

### Script Explanation

```bash
#!/bin/bash
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/secapr_env

# Define directories and input files
RAW_DIR="raw"
CLEAN_DIR="clean"
ASSEMBLY_DIR="assemblies"
BLAST_DIR="blast_results"
REF_FILE="B73_target_sequences.fasta"  
FINAL_HITS="final_combined_hits.fasta"
SAMPLE_FILE="sample_annotation.tab"  

# Create output directories if they don't exist
mkdir -p ${CLEAN_DIR}
mkdir -p ${ASSEMBLY_DIR}
mkdir -p ${BLAST_DIR}
mkdir -p logs

echo "Submitting BLAST jobs..."
while IFS=$'\t' read -r SAMPLE_ID FILENAME; do
    ASSEMBLY="${ASSEMBLY_DIR}/${SAMPLE_ID}/scaffolds.fasta"
    OUT_FILE="${BLAST_DIR}/${SAMPLE_ID}_blast_results.txt"
    BEST_HITS="${BLAST_DIR}/${SAMPLE_ID}_best_hits.fasta"
    
    # Submit BLAST job with dependency on SPAdes completion
    bsub -J "blast_${SAMPLE_ID}" \
         -o "logs/blast_${SAMPLE_ID}.out" \
         -e "logs/blast_${SAMPLE_ID}.err" \
         -n 4 -R "rusage[mem=8G]" \
         "# Wait until assembly file exists
          while [ ! -f ${ASSEMBLY} ]; do
              sleep 30
          done
          
          # Create BLAST database
          makeblastdb -in ${ASSEMBLY} -dbtype nucl -out ${ASSEMBLY_DIR}/${SAMPLE_ID}/db  -parse_seqids && \
          
          # Run BLAST and keep only best hit for each query
          blastn -task megablast \
                -query ${REF_FILE} -db ${ASSEMBLY_DIR}/${SAMPLE_ID}/db \
                -outfmt 6 \
                -num_alignments 1 -max_hsps 1 \
                -num_threads 4 > ${OUT_FILE} && \
          
          # Extract sequences for best hits
          # remove B73 from file name
          cut -f1,2 ${OUT_FILE} | sed s'/_B73//' | \
          while read QUERY SUBJECT; do \
              echo -e \">\${QUERY}_${SAMPLE_ID}\\n\$(blastdbcmd -db ${ASSEMBLY_DIR}/${SAMPLE_ID}/db -entry \${SUBJECT} -outfmt '%s')\"; \
          done > ${BEST_HITS}"
    
    echo "Submitted BLAST job for ${SAMPLE_ID}"
done < <(grep "_R1_" ${SAMPLE_FILE})


# Convert fasta files to one line sequence files
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf "%s",$0}}}' ../*_target_sequences.fasta  > target_sequences.fasta
# Make single file for all genes, all taxa
cat target_sequences.fasta   blast_results/*_best_hits.fasta | perl -pe 's/>/\n>/g' > blast_results/all_hits.fas 
# Group the sequences by gene
cut -f1 blast_results/*_blast_results.txt| sed 's/_.*//g'| sort |uniq | while IFS=$'\t' read -r GENE; do grep -A 1 "${GENE}" blast_results/all_hits.fas | grep -v -- "^--$" > blast_results/${GENE}.fas; done
```

### Running the Script

Ensure Step 2 has completed, then run:

```bash
chmod +x ./batch_scripts/q_find_target_contigs.sh
./batch_scripts/q_find_target_contigs.sh
```

### What Happens?

The script has two main parts:

1. **BLAST jobs for each sample:**
   - Creates a BLAST database from each sample's assembly
   - Runs BLAST to find the best match between reference sequences and the assembly
   - Extracts the matching sequences from the assembly
   - Outputs results to `blast_results/SAMPLE_ID_best_hits.fasta`

2. **Sequence organization by gene:**
   - Reformats reference sequences
   - Combines all hits and reference sequences
   - Organizes sequences by gene, creating separate files for each gene in `blast_results/GENE.fas`

### Key BLAST Parameters

- `-task megablast`: Uses megablast algorithm for highly similar sequences
- `-outfmt 6`: Outputs in tabular format
- `-num_alignments 1 -max_hsps 1`: Returns only the single best hit for each query

## Step 4: Multiple Sequence Alignment (q_align_sequences.sh)

This script creates multiple sequence alignments for each gene using MAFFT.

### Script Explanation

```bash
#!/bin/bash
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/secapr_env

# Define directories and input files
RAW_DIR="raw"
CLEAN_DIR="clean"
ASSEMBLY_DIR="assemblies"
BLAST_DIR="blast_results"
ALIGNMENT_DIR="alignments"  # New directory for alignments
REF_FILE="B73_target_sequences.fasta"  
FINAL_HITS="final_combined_hits.fasta"
SAMPLE_FILE="sample_annotation.tab"  

# Create output directories if they don't exist
mkdir -p ${CLEAN_DIR}
mkdir -p ${ASSEMBLY_DIR}
mkdir -p ${BLAST_DIR}
mkdir -p logs


# Add MAFFT alignment step by gene group
echo "Submitting MAFFT alignment jobs..."
cut -f1 ${BLAST_DIR}/*_blast_results.txt |sed 's/_.*//g'| sort | uniq | while IFS=$'\t' read -r GENE; do
    INPUT_FILE="${BLAST_DIR}/${GENE}.fas"
    OUTPUT_FILE="${ALIGNMENT_DIR}/${GENE}_aligned.fasta"
    
    # Create directory for the gene if it doesn't exist
    mkdir -p "${ALIGNMENT_DIR}/$(dirname ${GENE})"
    
    # Submit MAFFT job
    bsub -J "mafft_${GENE}" \
         -o "logs/mafft_${GENE}.out" \
         -e "logs/mafft_${GENE}.err" \
         -n 2 -R "rusage[mem=4G]" \
         -W 60 \
         "mafft --reorder --adjustdirection --nuc --auto  ${INPUT_FILE} | sed 's/^>_R_/>/' > ${OUTPUT_FILE}"
    echo "Submitted MAFFT alignment job for ${GENE}"
done
```

### Running the Script

Ensure Step 3 has completed, then run:

```bash
chmod +x ./batch_scripts/q_align_sequences.sh
./batch_scripts/q_align_sequences.sh
```

### What Happens?

- The script identifies all unique genes from the BLAST results
- For each gene, it submits a job to create a multiple sequence alignment using MAFFT
- The alignments are saved in the `alignments/` directory
- Each alignment contains all samples' sequences for a specific gene

### MAFFT Parameters

- `--reorder`: Outputs sequences in alignment order rather than input order
- `--adjustdirection`: Adjusts sequence direction to minimize mismatches
- `--nuc`: Specifies nucleotide sequences
- `--auto`: Automatically selects appropriate alignment strategy based on data

## Analysis and Visualization

After running the pipeline, you'll have:

1. Cleaned reads in `clean/`
2. De novo assemblies in `assemblies/`
3. BLAST results in `blast_results/`
4. Multiple sequence alignments in `alignments/`

These alignments can be used for:
- Phylogenetic analysis
- SNP identification
- Evolutionary studies
- Primer design

You can visualize the alignments using tools like:
- Jalview
- AliView
- MEGA
- R packages like `ape` or `Biostrings`

## Troubleshooting

### Common Issues

1. **Missing Dependencies**
   - Ensure your conda environment contains all necessary tools
   - Check for error messages about missing commands

2. **Job Failures**
   - Check the log files in the `logs/` directory
   - Look for resource constraints (time/memory)

3. **No BLAST Hits**
   - Verify that your reference sequences are appropriate for your samples
   - Check if assembly quality is sufficient

4. **LSF-Specific Issues**
   - If jobs are pending for a long time, check queue status with `bjobs`
   - Adjust resource requests if needed

### Log Files

Always check log files for errors:
```
cat logs/fastp_SAMPLE_ID.err
cat logs/spades_SAMPLE_ID.err
cat logs/blast_SAMPLE_ID.err
cat logs/mafft_GENE.err
```

## Customization Options

### Adjusting Resource Requests

If you need to adjust computing resources:

1. For read cleaning (Step 1):
   ```bash
   -n 4 -R "rusage[mem=4G]" -W 60
   ```

2. For assembly (Step 2):
   ```bash
   -n 8 -R "rusage[mem=20G]" -W 48:00
   ```

3. For BLAST (Step 3):
   ```bash
   -n 4 -R "rusage[mem=8G]"
   ```

4. For alignment (Step 4):
   ```bash
   -n 2 -R "rusage[mem=4G]"
   ```

### Quality Parameters

To adjust quality filtering parameters, modify the fastp command:
```bash
fastp -w 4 -i ${R1} -I ${R2} -o ${OUT_R1} -O ${OUT_R2} \
      --dedup \
      --qualified_quality_phred 30 \  # Adjust quality threshold
      --length_required 50 \          # Adjust minimum length
      ...
```

## Conclusion

This pipeline provides a streamlined approach for processing Illumina sequencing data, with each step building on the previous one. The modular design allows for troubleshooting at each stage and makes it easy to customize for different projects.

For maize genetics studies, this pipeline can help identify gene variants across different samples, potentially revealing important genetic differences related to traits of interest.