#Integration of h5 objs using seurat

setwd("D:/baylor/DRG atlas/integration/jungetal2023")

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)

source("utils/custom_seurat_functions.R")

GSM6069069 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069069_SAM24385747_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069070 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069070_SAM24385748_filtered_feature_bc_matrix.h5",use.names = T) 
GSM6069071 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069071_SAM24395829_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069072 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069072_SAM24395830_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069073 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069073_SAM24395831_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069074 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069074_SAM24395832_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069075 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069075_SAM24395833_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069076 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069076_SAM24395834_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069077 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069077_SAM24395835_filtered_feature_bc_matrix.h5",use.names = T)
GSM6069078 <- Read10X_h5("D:/baylor/DRG atlas/integration/jungetal2023/GSM6069078_SAM24395836_filtered_feature_bc_matrix.h5",use.names = T)

GSM6069069drg   <- CreateSeuratObject(GSM6069069, project = "GSM6069069")
GSM6069070drg   <- CreateSeuratObject(GSM6069070, project = "GSM6069070")
GSM6069071drg   <- CreateSeuratObject(GSM6069071, project = "GSM6069071")
GSM6069072drg   <- CreateSeuratObject(GSM6069072, project = "GSM6069072")
GSM6069073drg   <- CreateSeuratObject(GSM6069073, project = "GSM6069073")
GSM6069074drg   <- CreateSeuratObject(GSM6069074, project = "GSM6069074")
GSM6069075drg   <- CreateSeuratObject(GSM6069075, project = "GSM6069075")
GSM6069076drg   <- CreateSeuratObject(GSM6069076, project = "GSM6069076")
GSM6069077drg   <- CreateSeuratObject(GSM6069077, project = "GSM6069077")
GSM6069078drg   <- CreateSeuratObject(GSM6069078, project = "GSM6069078")

rm(GSM6069069)
rm(GSM6069070)
rm(GSM6069071)
rm(GSM6069072)
rm(GSM6069073)
rm(GSM6069074)
rm(GSM6069075)
rm(GSM6069076)
rm(GSM6069077)
rm(GSM6069078)


GSM6069069drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069069drg, pattern = "^MT-")
GSM6069069drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069069drg, pattern = "^Mt-")
VlnPlot(GSM6069069drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

GSM6069070drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069070drg, pattern = "^MT-")
GSM6069070drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069070drg, pattern = "^Mt-")
VlnPlot(GSM6069070drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

GSM6069071drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069071drg, pattern = "^MT-")
GSM6069071drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069071drg, pattern = "^Mt-")
VlnPlot(GSM6069071drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

GSM6069072drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069072drg, pattern = "^MT-")
GSM6069072drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069072drg, pattern = "^Mt-")
VlnPlot(GSM6069072drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

GSM6069073drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069073drg, pattern = "^MT-")
GSM6069073drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069073drg, pattern = "^Mt-")
VlnPlot(GSM6069073drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

GSM6069074drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069074drg, pattern = "^MT-")
GSM6069074drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069074drg, pattern = "^Mt-")
VlnPlot(GSM6069074drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)


GSM6069075drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069075drg, pattern = "^MT-")
GSM6069075drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069075drg, pattern = "^Mt-")
VlnPlot(GSM6069075drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)


GSM6069076drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069076drg, pattern = "^MT-")
GSM6069076drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069076drg, pattern = "^Mt-")
VlnPlot(GSM6069076drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

GSM6069077drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069077drg, pattern = "^MT-")
GSM6069077drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069077drg, pattern = "^Mt-")
VlnPlot(GSM6069077drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

GSM6069078drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069078drg, pattern = "^MT-")
GSM6069078drg[["percent.mt"]]  <- PercentageFeatureSet(GSM6069078drg, pattern = "^Mt-")
VlnPlot(GSM6069078drg, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)


drg_list <- list()
drg_list[["GSM6069069drg"]] <- GSM6069069drg
drg_list[["GSM6069070drg"]] <- GSM6069070drg
drg_list[["GSM6069071drg"]] <- GSM6069071drg
drg_list[["GSM6069072drg"]] <- GSM6069072drg
drg_list[["GSM6069073drg"]] <- GSM6069073drg
drg_list[["GSM6069074drg"]] <- GSM6069074drg
drg_list[["GSM6069075drg"]] <- GSM6069075drg
drg_list[["GSM6069076drg"]] <- GSM6069076drg
drg_list[["GSM6069077drg"]] <- GSM6069077drg
drg_list[["GSM6069078drg"]] <- GSM6069078drg

for (i in 1:length(drg_list)) {
  drg_list[[i]] <- NormalizeData(drg_list[[i]], verbose = F)
  drg_list[[i]] <- FindVariableFeatures(drg_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = F)
}

drg_anchors <- FindIntegrationAnchors(object.list = drg_list, dims = 1:30)
drg_seurat  <- IntegrateData(anchorset = drg_anchors, dims = 1:30)

rm(drg_list)
rm(drg_anchors)

DefaultAssay(drg_seurat) <- "RNA"

DefaultAssay(drg_seurat) <- "integrated"
drg_seurat <- ScaleData(drg_seurat, verbose = F)
drg_seurat <- RunPCA(drg_seurat, npcs = 30, verbose = F)
drg_seurat <- RunUMAP(drg_seurat, reduction = "pca", dims = 1:30, verbose = F)

DimPlot(drg_seurat, reduction = "umap") + plot_annotation(title = "jungetal2023 INTEGRATED")

DimPlot(drg_seurat, reduction = "umap", split.by = "orig.ident") + NoLegend()

drg_seurat <- FindNeighbors(drg_seurat, dims = 1:30, k.param = 10, verbose = F)
drg_seurat <- FindClusters(drg_seurat, verbose = F)
DimPlot(drg_seurat,label = T) + NoLegend()

count_table <- table(drg_seurat@meta.data$seurat_clusters, drg_seurat@meta.data$orig.ident)
count_table

plot_integrated_clusters(drg_seurat) 

saveRDS(drg_seurat, "jungetal2023.rds")
