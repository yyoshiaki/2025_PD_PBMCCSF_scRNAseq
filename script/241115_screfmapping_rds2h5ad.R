library(zellkonverter)
library(Seurat)

setwd('/home/yy693/pi_hafler/ASAP')

for (dir.input in  list.dirs('output/cellranger', recursive = FALSE)) {
    print(paste('subdir: ', dir.input))
    # # take the sample name from the command line
    # subdir <- commandArgs(trailingOnly = TRUE)[1]
  
  # samples <- c('YPD001BLD')
  # dir.input <- paste0("output/cellranger/", subdir)
  # get the list of samples in the directory
  # samples <- list.files(dir.input)
  
  sample_dirs <- list.dirs(path = dir.input, recursive = FALSE, full.names = TRUE)
  
  valid_samples <- sapply(sample_dirs, function(dir) {
    file.exists(file.path(dir, "outs", "filtered_feature_bc_matrix.h5"))
  })
  
  valid_sample_dirs <- sample_dirs[valid_samples]
  
  samples <- basename(valid_sample_dirs)
  
  
  ## ----parameter setting (change here)---------------------------------------------------------------
  for (project.name in samples) {
    print(paste('start analysis: ', project.name))
    prefix <- paste0('output/CD4T_screfmapping/', basename(dir.input), '/', project.name, "/", project.name)
    
    if (file.exists(paste0(prefix, "_CD4T_AssayData.rds")) && 
        !file.exists(paste0(prefix, "_CD4T_AssayData.h5ad"))) {
      
      print(paste('Making H5AD: ', project.name))
      # load extracted CD4T
      query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))
      query_obj <- CreateSeuratObject(counts = query_obj,
                                      project = project.name,
                                      assay = "RNA",
                                      min.cells = 3,
                                      min.features = 200)
      writeH5AD(as.SingleCellExperiment(query_obj), 
                paste0(prefix, "_CD4T_AssayData.h5ad"))
      print(paste('H5AD saved: ', project.name))
    }
  }
}
