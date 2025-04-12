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
`fastp`  from https://github.com/OpenGene/fastp

```
conda install -c bioconda fastp
```

`circos` has `perl::GD` for plotting of blast results.

```
conda create -p /share/maize/user/env/circos  bioconda::circos bioconda::blast
conda activate /share/maize/user/env/circos
```

`blast-imager.pl`
https://github.com/vinuesa/TIB-filoinfo/tree/master

