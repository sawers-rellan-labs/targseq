#!/bin/bash
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/secapr_env

# Define directories and input files
RAW_DIR="raw"
CLEAN_DIR="clean"
ASSEMBLY_DIR="assemblies"
BLAST_DIR="blast_results"
REF_FILE="reference_sequences.fasta"  # Your reference file with REF01-REF16
FINAL_HITS="final_combined_hits.fasta"
SAMPLE_FILE="sample_annotation.tab"  # Tab-separated file with sample info

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