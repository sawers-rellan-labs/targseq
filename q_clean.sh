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
             "fastp -i ${R1} -I ${R2} -o ${OUT_R1} -O ${OUT_R2} \
                    --detect_adapter_for_pe \
                    --correction \
                    --cut_front \
                    --cut_tail \
                    --cut_window_size 4 \
                    --cut_mean_quality 20 \
                    --qualified_quality_phred 20 \
                    --length_required 50 \
                    --html ${CLEAN_DIR}/${SAMPLE_ID}_fastp.html \
                    --json ${CLEAN_DIR}/${SAMPLE_ID}_fastp.json"

        echo "Submitted fastp job for ${SAMPLE_ID} (${BASE_FILENAME})"
    else
        echo "Warning: Input files for ${SAMPLE_ID} not found. Skipping."
    fi
# Only process lines that have _R1_ in them to avoid duplicates (since file contains both R1 and R2)
done < <(grep "_R1_" ${SAMPLE_FILE})