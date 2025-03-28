## Set up
We installed `secaper` from:
https://github.com/AntonelliLab/seqcap_processor

See `secaper_install.md`

#### Login to hazel
```
ssh user@login.hcp.ncsu.edu
```

#### Request an interactive session
```
bsub -q sara -Is -n 8 -R "rusage[mem=128]" -W 8:00 bash
```

#### Activate `secaper` conda environment
```
conda activate /usr/local/usrapps/maize/user/secapr_env
```

#### Directory structure
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
#### Create a subset of data for testing
```
mkdir testin
ls -1  data/*fastq.gz | tail -n 8 | xargs -I{} cp {} testin
```

#### Execute code from your user folder
```
cd user
mkdir testout
```

## Tutorial
We are following the tutorial from with our own data

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
