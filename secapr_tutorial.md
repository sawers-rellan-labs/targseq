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
tree ../targseq
```

#### Create a subset of data for testing
```
mkdir testin
ls -1  data/*fastq.gz | tail -n 8 | xargs -I{} cp {} testin
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
```
secapr clean_reads --input testin  --output testout
```

### Reference assembly
