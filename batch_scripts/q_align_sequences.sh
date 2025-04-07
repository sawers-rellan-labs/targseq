#!/bin/bash
source ~/.bashrc
conda activate /share/maize/frodrig4/conda/env/secapr_env

# Define directories and input files
RAW_DIR="raw"
CLEAN_DIR="clean"
ASSEMBLY_DIR="assemblies"
BLAST_DIR="blast_results"
ALIGNMENT_DIR="alignments"  # New directory for alignments
REF_FILE="reference_sequences.fasta"  # Your reference file with REF01-REF16
FINAL_HITS="final_combined_hits.fasta"
SAMPLE_FILE="sample_annotation.tab"  # Tab-separated file with sample info

# Create output directories if they don't exist
mkdir -p ${CLEAN_DIR}
mkdir -p ${ASSEMBLY_DIR}
mkdir -p ${ALIGNMENT_DIR}
mkdir -p ${BLAST_DIR}
mkdir -p "${BLAST_DIR}/by_gene"
mkdir -p logs


# Group the sequences by gene
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf "%s",$0}}}' ../reference_sequences.fasta  > reference_sequences.fasta
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf "%s",$0}}}' ../TIL18_reference_sequences.fasta  > TIL18_reference_sequences.fasta
cat reference_sequences.fasta TIL18_reference_sequences.fasta  blast_results/*_best_hits.fasta | perl -pe 's/>/\n>/g' > blast_results/all_hits.fas 

cut -f1 blast_results/*_blast_results.txt| sort |uniq | sed 's/_.*//'| while IFS=$'\t' read -r GENE; do grep -A 1 "${GENE}" blast_results/all_hits.fas | grep -v -- "^--$" > blast_results/by_gene/${GENE}.fas; done

# Add MAFFT alignment step by gene group
echo "Submitting MAFFT alignment jobs..."
cut -f1 ${BLAST_DIR}/*_blast_results.txt | sort | uniq | sed 's/_.*//'| while IFS=$'\t' read -r GENE; do
    INPUT_FILE="${BLAST_DIR}/by_gene/${GENE}.fas"
    OUTPUT_FILE="${ALIGNMENT_DIR}/${GENE}_aligned.fasta"
    
    # Create directory for the gene if it doesn't exist
    mkdir -p "${ALIGNMENT_DIR}/$(dirname ${GENE})"
    
    # Submit MAFFT job
    bsub -J "mafft_${GENE}" \
         -o "logs/mafft_${GENE}.out" \
         -e "logs/mafft_${GENE}.err" \
         -n 2 -R "rusage[mem=4G]" \
         "mafft --reorder --adjustdirection --nuc --auto  ${INPUT_FILE} | sed 's/^>_R_/>/' > ${OUTPUT_FILE}"
    echo "Submitted MAFFT alignment job for ${GENE}"
done

