setwd("D:/baylor/TRIGEMINAL/integration/b6tgnu_neuron10x")

# To load the libraries

library(Seurat)
library(SeuratData)
library(patchwork)

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

# To Identify conserved cell type markers

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

# To generate the feature plot to a given gene list

FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                          "CCL2", "PPBP"), min.cutoff = "q9")


#################################################################################################################
#                                  Merging More Than Two Seurat Objects                                         #
#################################################################################################################

# First we need to load the .h5 files

test1 <- Read10X_h5("D:/baylor/TRIGEMINAL/B6_TG_Nu/analysis/input/filtered_feature_bc_matrix.h5")
test2 <- Read10X_h5("D:/baylor/TRIGEMINAL/cellqcmar5/GSM5912707/GSM5912707.h5")
test3 <- Read10X_h5("D:/baylor/TRIGEMINAL/B6_TG_Nu/integration/b6tgnu_neuron10x/GSM5912706.h5")
test4 <- Read10X_h5("D:/baylor/TRIGEMINAL/B6_TG_Nu/integration/b6tgnu_neuron10x/GSM5967143.h5")


test.big <- merge(test1, y = c(test2, test3, test4), add.cell.ids = c("6K", "10K", "8K", "7K"), project = "testintegrateobj")

test.big











