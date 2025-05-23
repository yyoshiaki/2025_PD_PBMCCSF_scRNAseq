## prep sif file

3/22/2024

```
cd /home/yy693/pi_hafler/Tetramer_scRNAseq
apptainer build screfmapping-0.0.1.sif docker://yyasumizu/screfmapping:0.0.1
```
modified to `k.anchor = 3` in '/home/yy693/pi_hafler/Tetramer_scRNAseqscrefmapping/ref_mapping_seuratobj.R' as in 
https://github.com/yyoshiaki/screfmapping/tree/main?tab=readme-ov-file#the-number-of-neighbors-k-to-use-when-finding-anchors

```
docker run --rm -it -v ${PWD}:/home/rstudio/autoimmune_10x  yyasumizu/screfmapping:0.0.1 Rscript script/screfmapping.R

apptainer exec screfmapping-0.0.1.sif Rscript script/screfmapping.R
```