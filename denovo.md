```
ls -1 ../testin/ | \
   perl -pe 's/(^.*(S\d+))_L001/$2,$1\_L001/' | \
   sed 's/\.fastq\.gz//g' \
   > sample_annotation.csv
```   

```
mkdir raw
ls -1 ../data | grep -P "_S20_|_S27_|_S29_" | xargs -I{} cp ../data/{} raw

ls -1 raw/ |    perl -pe 's/(^.*(S\d+))_L001/$2\t$1\_L001/'| sort -n -k1.2 > sample_annotation.tab
```