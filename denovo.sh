# for 
fastp -i ../testin/Zx0540_P3_P5_P1111_S27_L001_R1_001.fastq.gz \
      -o Zx0540_P3_P5_P1111_S27_L001_R1_clean.fastq.gz \
      -I ../testin/Zx0540_P3_P5_P1111_S27_L001_R2_001.fastq.gz \
      -O Zx0540_P3_P5_P1111_S27_L001_R2_clean.fastq.gz 

# probably I can make this with bsub on the sample name

# using fastp for cleaning all the reads
# for R1 in testin/*_R1_001.fastq.gz; do R2=$(echo $R1| sed 's/_R1_/_R2_/'); name=$(echo $R1|sed 's/_R1_001.fastq.gz//'|sed 's/Raw_Reads\///'); echo "fastp -i $R1 -I $R2 -o clean/$name.clean.R1.fastq.gz -O clean/$name.clean.R2.fastq.gz -w 4 --detect_adapter_for_pe -j log/clean/$name.json -h log/clean/$name.html &> log/clean/$name.log.txt";done > fastp.commands

# fastqc -t 48 -o fastqc_output/ *.gz
# cd fastqc_output
# multiqc .

# this did not work
# secapr assemble_reads  --kmer 35 --input clean --output contigs 


spades.py -t 8 -k 35\
          --isolate \
          -1 clean/S27/Zx0540_P3_P5_P1111_S27_L001_R1_clean.fastq.gz \
          -2 clean/S27/Zx0540_P3_P5_P1111_S27_L001_R2_clean.fastq.gz \
          -o contigs

# │   ├── configs
# │   │   ├── careful_mda_mode.info
# │   │   ├── careful_mode.info
# │   │   ├── config.info
# │   │   ├── construction.info
# │   │   ├── detail_info_printer.info
# │   │   ├── distance_estimation.info
# │   │   ├── hmm_mode.info
# │   │   ├── isolate_mode.info
# │   │   ├── large_genome_mode.info
# │   │   ├── mda_mode.info
# │   │   ├── meta_mode.info
# │   │   ├── metaplasmid_mode.info
# │   │   ├── metaviral_mode.info
# │   │   ├── pe_params.info
# │   │   ├── plasmid_mode.info
# │   │   ├── rna_mode.info
# │   │   ├── rnaviral_mode.info
# │   │   ├── sewage_mode.info
# │   │   ├── simplification.info
# │   │   └── toy.info


contig_stats.pl contigs/contigs.fasta

