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

#### Directory structure
```
cd /rsstu/users/r/rrellan/BZea/targseq/
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
ls  -aFlhSr data| head -n 20
```

select some megabase size files for the test

```
mkdir testin
ls -1 data | grep -P "_S20_|_S27_|_S29_" | xargs -I{} cp data/{} testin
```




```
tree ../targseq
```

```
.
|--|data
|--|testin
|--|user
    |--|testout
    
```


#### Execute code from your user folder
```
cd user
mkdir testout
```

## Tutorial
We are following the [secapr tutorial](https://htmlpreview.github.io/?https://github.com/AntonelliLab/seqcap_processor/blob/master/docs/documentation/tutorial.html)
with our own data

### *De novo* assembly 

#### Quality check
```
secapr quality_check --input  ../testin/ --output testout/qc
```


#### Clean reads

Make sample annotation table
```
ls -1 ../testin/ | \
   perl -pe 's/^(.*?)_L001/$1,$1\_L001/' | \
   sed 's/\.fastq\.gz//g' \
   > sample_annotation.csv
```

Run clean command 
```
secapr clean_reads \
       --input ../testin \
       --read_min 190000 \
       --sample_annotation_file sample_annotation.csv \
       --output testout/clean
```

### Reference assembly
