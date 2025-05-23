library(Seurat)
library(MuDataSeurat)

obj <- readRDS('./data/PDnuclei.rds')
obj
WriteH5AD(obj, "./data/PDnuclei.h5ad", assay = "RNA")
write.csv(obj@reductions[["umap"]]@cell.embeddings, file = './data/PDnuclei.umap.csv')