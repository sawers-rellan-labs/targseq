## Set up
We installed `secapr` from:
https://github.com/AntonelliLab/seqcap_processor

See `secapr_installation.md`

#### Login to hazel
```
ssh user@login.hpc.ncsu.edu
```

#### Request an interactive session
```
bsub -q sara -Is -n 8 -R "rusage[mem=128]" -W 8:00 bash
```

#### Activate `secapr` conda environment
```
conda activate /usr/local/usrapps/maize/user/secapr_env
```


#### Execute code from your user folder
```
cd /rsstu/users/r/rrellan/BZea/targseq/
cd user
mkdir testin
mkdir testout
```

#### Directory structure
```
tree ../targseq
```

#### Create a subset of data for testing

```
# selecting first four samples for test
# ls -1  data/*fastq.gz | tail -n 8 | xargs -I{} cp {} testin
# Didn't work files too large to test
```

check for small files
```
ls  -aFlhSr ../data| head -n 20
```

select some megabase size files for the test

```
mkdir testin
ls -1 ../data | grep -P "_S20_|_S27_|_S29_" | xargs -I{} cp ../data/{} testin
# ls -1 ../data | grep -P "_S27_|_S1_|_S9_" | xargs -I{} cp ../data/{} raw
```



```
tree ./
```

```

|--|user
    |--|testout
    |--|testin
```


## Tutorial
We are following the [secapr tutorial](https://htmlpreview.github.io/?https://github.com/AntonelliLab/seqcap_processor/blob/master/docs/documentation/tutorial.html)
with our own data

### *De novo* assembly 

#### Quality check
```
secapr quality_check --input  testin/ --output testout/qc
```

#### Clean reads
Make sample annotation table is a table relateing sample names to the sastqc files.
```
Zx0540_P3_P5_P1111_S27,Zx0540_P3_P5_P1111_S27_L001_R1_001.fastq
Zx0540_P3_P5_P1111_S27,Zx0540_P3_P5_P1111_S27_L001_R2_001.fastq
Zx0550_P4_P3_P5311_S20,Zx0550_P4_P3_P5311_S20_L001_R1_001.fastq
Zx0550_P4_P3_P5311_S20,Zx0550_P4_P3_P5311_S20_L001_R2_001.fastq
Zx0580_P2_P5_P2411_S29,Zx0580_P2_P5_P2411_S29_L001_R1_001.fastq
Zx0580_P2_P5_P2411_S29,Zx0580_P2_P5_P2411_S29_L001_R2_001.fastq
```
You can make ut ias a spreadsheet and save it as csv.
But we'll run a shell command using perl and regex to generate the annotation file from the folser contents.
```
ls -1 testin/ | \
   perl -pe 's/^(.*?)_L001/$1,$1\_L001/' | \
   sed 's/\.fastq\.gz//g' \
   > sample_annotation.csv
```


The originial command in the tutorial does nopt work because it is an old version.

```
# secapr clean_reads --input pipeline_exercise/fastq_raw/ --config pipeline_exercise/adapter_info.txt --output pipeline_exercise/cleaned_trimmed_reads --index singlepipeline_exercise/fastqc_results/raw
```
That version needed `adapter_info.txt`, but tyhe current version without it.

This new command is supossed to work

Run clean command 
```
secapr clean_reads \
       --input testin \
       --read_min 190000 \
       --sample_annotation_file sample_annotation.csv \
       --output clean
```


it is not working see

https://github.com/AntonelliLab/seqcap_processor/issues/41


### Reference assembly

```
secapr reference_assembly --reads pipeline_exercise/cleaned_trimmed_reads --reference_type alignment-consensus --reference pipeline_exercise/alignments/contig_alignments --output pipeline_exercise/mapped_reads --min_coverage 4
```

```
secapr locus_selection --input pipeline_exercise/mapped_reads --output pipeline_exercise/selected_loci --n 100
```

```
secapr align_sequences --sequences pipeline_exercise/mapped_reads/joined_unphased_fastas.fasta --outdir pipeline_exercise/alignments/bam_consensus_alignments --no_trim
```

```
secapr phase_alleles --input pipeline_exercise/selected_loci --output pipeline_exercise/allele_sequences_selected_loci --min_coverage 3
```

```
secapr align_sequences --sequences pipeline_exercise/allele_sequences_selected_loci/joined_allele_fastas.fasta --outdir pipeline_exercise/alignments/selected_loci_allele_alignments --no_trim
```




