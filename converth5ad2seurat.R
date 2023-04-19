#To convert .h5ad to h5Seurat and to h5 obj
setwd("D:/baylor/TRIGEMINAL/B6_TG_Nu/integration/b6tgnu_neuron10x")
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(reticulate)
library(anndata)
# To load the .h5ad file
objintegrated <- anndata::read_h5ad("D:/baylor/TRIGEMINAL/B6_TG_Nu/integration/b6tgnu_neuron10x/objintegrated.h5ad")
Convert("objintegrated.h5ad", dest = "h5seurat", overwrite = TRUE)
objintegrated <- LoadH5Seurat("objintegrated.h5seurat")
objintegrated
