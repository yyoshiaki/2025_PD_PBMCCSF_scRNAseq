# 2025 PD PB & CSF scRNA-seq Integrated Analysis

## Overview
This repository contains notebooks and scripts for preprocessing, clustering, sub-clustering, and integrative analysis of single-cell RNA sequencing (scRNA-seq) data derived from peripheral blood (PB) and cerebrospinal fluid (CSF) of Parkinson’s disease (PD) patients.

## Directory Structure

### Preprocess

- notebook/
    - 241006_scvi_BLD.py  
    - 241006_scvi_CSF.py  
    - 241010_scvi_BLD_subcluster_B.py  
    - 241010_scvi_BLD_subcluster_myelo.py  
    - 241010_scvi_BLD_subcluster_NK.py  
    - 241010_scvi_BLD_subcluster_T.py  
    - 241026_scvi_CSF_subcluster_B.py  
    - 241026_scvi_CSF_subcluster_myelo.py  
    - 241026_scvi_CSF_subcluster_T.py  
    - 241026_scvi_CSF_subcluster_TNK.py  
    - 241021_BLD_integrate.ipynb  
    - 241027_CSF_integrate.ipynb  
### Downstream analysis
- notebook/
    - 250407_PBMC_downstream.ipynb
    - 250407_CSF_downstream.ipynb
    - 241114_CSF_myelo_downstream.ipynb  
    - 250407_scCODA_BLD.ipynb
    - 250407_scCODA_CSF.ipynb
### Comparison with MS
- notebook/
    - 250417_MS_scanvi.py
    - 250417_CSFMac_MS_downstream.ipynb

### Misc
- notebook/
    - 241113_CSF_scDRS.ipynb  : scDRS
    - 241127_CSF_cpdb.ipynb  : cellphone DB
    - 250501_mouseBAM_Dennis.ipynb : mouse dura

### VDJ analysis
- notebook/
    - 241212_BLD_VDJ.ipynb
    - 241212_CSF_VDJ.ipynb
    - 241219_startrac.ipynb
- script/  

## dependencies

- scanpy==1.10.2
- scirpy==0.20.0
- scvi-tols==1.1.6

## License
This project is licensed under the [MIT License](LICENSE).  

## Citation

Zhang, Le, Yoshiaki Yasumizu, Marion Elizabeth Deerhake, Jeonghyeon Moon, Nicholas Buitrago-Pocasangre, Anthony Russo, Haowei Wang, et al. 2025. “Inflamed Microglia like Macrophages in the Central Nervous System of Prodromal Parkinson′s Disease.” bioRxiv. https://doi.org/10.1101/2025.05.16.654530.

## Acknowledgements

This research was funded in part by Aligning Science Across Parkinson’s [ASAP-000529] through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).
