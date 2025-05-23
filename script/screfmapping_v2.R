## base: example.R

## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)
library(Matrix)

# setwd('/home/rstudio/autoimmune_10x')
setwd('/home/yy693/pi_hafler/ASAP')
source('/home/yy693/pi_hafler/ASAP/screfmapping/ref_mapping_seuratobj.R')
source('/home/yy693/pi_hafler/ASAP/screfmapping/utils_seurat.R')

options(Seurat.object.assay.calcn = TRUE) 

# subdir <- "Control"
# take the sample name from the command line
subdir <- commandArgs(trailingOnly = TRUE)[1]

# samples <- c('YPD001BLD')
dir.input <- paste0("output/cellranger/", subdir)
# get the list of samples in the directory
# samples <- list.files(dir.input)

sample_dirs <- list.dirs(path = dir.input, recursive = TRUE, full.names = TRUE)

valid_samples <- sapply(sample_dirs, function(dir) {
  file.exists(file.path(dir, "outs", "filtered_feature_bc_matrix.h5"))
})

valid_sample_dirs <- sample_dirs[valid_samples]

samples <- basename(valid_sample_dirs)

## ----Load　reference---------------------------------------------------------------------------
# Reference for CD4T cells are included in Docker repo. No need to modify this section for CD4T analysis.
# Azimuth
# reference <- LoadReference(path = "./screfmapping/data/Azimuth/human_pbmc_v1.0.0")
reference <- LoadReference(path = "/screfmapping/data/Azimuth/human_pbmc_v1.0.0")
# Symphony
# load("./screfmapping/data/ref_Reference_Mapping_20220525.RData")
load("/screfmapping/data/ref_Reference_Mapping_20220525.RData")

file.copy(from = '/screfmapping/data/cache_symphony_sct.uwot', 
  to = '/home/yy693/pi_hafler/ASAP/cache_symphony_sct.uwot')

## ----parameter setting (change here)---------------------------------------------------------------
for (project.name in samples) {
    print(paste('start analysis: ', project.name))
    prefix <- paste0("./output/CD4T_screfmapping/", subdir, '/', project.name, "/", project.name)

    # skip if the ouput already exists
    if (file.exists(paste0(prefix, "_Reference_Mapping.csv"))) {
        print(paste('skip analysis, output already exists: ', project.name))
        next
    }

    dir.create(paste0("./output/CD4T_screfmapping/", subdir, '/', project.name), recursive = T, showWarnings = FALSE)

    ## ----reading data (change here)--------------------------------------------------------------------
    # Load the PBMC dataset
    pbmc.data <- Read10X(data.dir = paste0("output/cellranger/", subdir, '/', project.name, '/outs/filtered_feature_bc_matrix'))
    q <- CreateSeuratObject(counts = pbmc.data,
                            project = project.name,
                            assay = "RNA",
                            min.cells = 3,
                            min.features = 200)
    
    print(q)
    
    # if object q conteins 0 cells, skip the analysis
    if (dim(q)[2] < 200) {
        print(paste('skip analysis, less than 200 cells detected: ', project.name))
        next
    }

    ## ----extraction of CD4T-----------------------------------------------------------------------------
    extract_cells_seuratobj(q, reference, prefix)

    # load extracted CD4T
    query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))

    if (dim(query_obj)[2] < 10) {
        print(paste('skip analysis, less than 10 CD4T cells detected: ', project.name))
        next
    }

    query_obj <- CreateSeuratObject(counts = query_obj,
                            project = project.name,
                            assay = "RNA",
                            min.cells = 3,
                            min.features = 200)
    
    # writeH5AD(as.SingleCellExperiment(query_obj), paste0(prefix, "_CD4T_AssayData.h5ad"))

    ## ----run Symphony-----------------------------------------------------------------------------
    reference_mapping_seuratobj(ref, query_obj, prefix)

}
## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
