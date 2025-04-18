
# Targseq
Bzea target sequencing of HPC1 and nitrogen genes

# Table of contents

- [Targseq](#targseq)
  - [Part 1: Get Target sequences](#part-1-get-target-sequences)
  - [Part 2: De novo assembly of target sequencing data](#part-2-de-novo-assembly-of-target-sequencing-data)

## Part 1: Get Target sequences
```mermaid
flowchart TD
    A[Start: B73 gene targets] --> B[Download ref genomes]
    B --> C[Create BLAST databases]
    
    C --> D{OrthoMCL available?}
    D -->|Yes| E1[Extract orthologs from OrthoMCL]
    D -->|No| E2[BLAST against each genome]
    
    E1 --> F[Create Ortholog Table]
    E2 --> F
    
    F --> H[Extract canonical transcripts]
    H --> I[Add 2kb flanking regions]
    I --> J[Create taxa database file]
    J --> K[Extract target sequences with blastdbcmd]
    K --> L[Rename sequences with gene symbols]
    L --> M[Combine sequences into final FASTA files]
    
    subgraph "Input Files"
    A1[B73_gene_targets.tab]
    end
    
    subgraph "Reference Data"
    B1[Fasta Genomes]
    B2[GFF3 Annotations]
    end
    
    subgraph "Output Files"
    O1[GENOME_target_sequences.fasta Target sequences with flanking regions]
    end
    
    A1 --> A
    B1 --> B
    B2 --> B

    
    M --> O1
``` 

## Part 2: De novo assembly of target sequencing data

```mermaid
flowchart TD
    subgraph "Setup"
        A[Raw Illumina Reads] --> |sample_annotation.tab| B["Directory Setup
        - raw/
        - clean/ 
        - assemblies/
        - blast_results/
        - alignments/
        - logs/"]
        C["Reference Files
        - B73_target_sequences.fasta
        - TIL18_target_sequences.fasta"] --> B
    end
  
    subgraph "Step 1: Read Cleaning"
        B --> D["q_clean_reads.sh"]
        D --> |"fastp
        - dedup
        - quality phred ≥30
        - length ≥50bp"| E["Cleaned Reads
        sample_R1_001.clean.fastq.gz
        sample_R2_001.clean.fastq.gz"]
    end
    
    subgraph "Step 2: De Novo Assembly"
        E --> F["q_assemble_reads.sh"]
        F --> |"SPAdes
        - isolate mode
        - auto coverage cutoff
        - 8 threads, 20GB memory"| G["De Novo Assemblies
        scaffolds.fasta"]
    end
    
    subgraph "Step 3: BLAST Search"
        G --> H["q_find_target_contigs.sh"]
        C --> H
        H --> |"BLAST
        - makeblastdb
        - blastn megablast
        - extract best hits"| I["Target Contigs
        sample_best_hits.fasta
        gene.fas"]
    end
    
    subgraph "Step 4: Multiple Sequence Alignment"
        I --> J["q_align_sequences.sh"]
        J --> |"MAFFT
        - reorder
        - adjustdirection
        - auto algorithm"| K["Multiple Sequence Alignments
        gene_aligned.fasta"]
    end
```
## Part 3: Phylogeny estimation

```mermaid
flowchart TD
    subgraph "Setup"
        A["Input Files
        - hpc1_aligned.fasta (Alignment)
        - seqid_label.csv (Metadata)"] --> B["Environment Setup
        - BiocManager
        - Biostrings
        - tidyverse
        - ape
        - phangorn
        - phytools
        - ggtree"]
    end
    
    subgraph "Step 1: Data Processing"
        B --> C["Read Input Files
        - Read alignment (readDNAStringSet)
        - Read metadata (read.csv)"]
        C --> D["Create Name Mapping
        - Extract sequence IDs
        - Join with metadata
        - Filter to match sequences"]
        D --> E["Relabel Sequences
        - Replace technical IDs with
        meaningful taxonomic labels"]
        E --> F["Trim Alignment
        - Focus on region starting 
        from position 2226"]
    end
    
    subgraph "Step 2: Variant Analysis"
        F --> G["Extract Variant Information
        - A204T variant
        - I211V variant
        - Convert to REF/ALT format"]
        G --> H["Write Processed Alignment
        - Save relabeled & trimmed alignment"]
    end
    
    subgraph "Step 3: Phylogenetic Analysis"
        H --> I["Build Full Phylogenetic Tree
        - Convert to phyDat format
        - Calculate distance matrix (JC69)
        - Build UPGMA tree"]
        I --> J["Tree Manipulation
        - Rotate tree (B73 at top)
        - Save tree in Newick format"]
    end
    
    subgraph "Step 4: Full Tree Visualization"
        J --> K["Create Variant Heatmap Tree
        - Add variant data
        - Add ancestry information
        - Color tips by ancestry"]
        K --> L["Add Multiple Alignment
        - Combine tree with alignment view
        - Save as PDF and PNG"]
    end
    
    subgraph "Step 5: Subset Analysis"
        J --> M["Create Donor Subset
        - Filter for donor ancestry + B73
        - Extract subset alignment
        - Build new tree for subset"]
        M --> N["Visualize Donor Subset
        - Create variant heatmap
        - Add alignment visualization
        - Save as PDF and PNG"]
    end
    
    style A fill:#f9f9f9,stroke:#333,stroke-width:1px
    style F fill:#e6f3ff,stroke:#333,stroke-width:1px
    style G fill:#e6f3ff,stroke:#333,stroke-width:1px
    style J fill:#e6f3ff,stroke:#333,stroke-width:1px
    style L fill:#f9fff9,stroke:#333,stroke-width:1px
    style N fill:#f9fff9,stroke:#333,stroke-width:1px
```  