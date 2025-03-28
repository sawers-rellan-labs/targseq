### Install

We installed `secapr` from:
https://github.com/AntonelliLab/seqcap_processor


#### Login to hazel
```
ssh user@login.hpc.ncsu.edu
```

#### Installing with `conda`
This step should be run from the login node.
```
 conda create -p /share/maize/user/conda/env/secapr_env secapr
 conda activate /share/maize/user/conda/env/secapr_env
```

#### Additional packages
`fasp` 
https://github.com/OpenGene/fastp

```
 conda install -c bioconda fastp
```
 