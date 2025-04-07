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
          cut -f1,2 ${OUT_FILE} | sed s'/_B73//'| \
          while read QUERY SUBJECT; do \
              echo -e \">\${QUERY}_${SAMPLE_ID}\\n\$(blastdbcmd -db ${ASSEMBLY_DIR}/${SAMPLE_ID}/db -entry \${SUBJECT} -outfmt '%s')\"; \
          done > ${BEST_HITS}"
    
    echo "Submitted BLAST job for ${SAMPLE_ID}"
done < <(grep "_R1_" ${SAMPLE_FILE})
