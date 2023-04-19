setwd("D:/baylor/HumanAnteriorSegment/fetalHAS")
library(Seurat)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)
library(mclust)
# load data

obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- Read10X_h5(filename = "D:/baylor/TRIGEMINAL/obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined/analysis/input/filtered_feature_bc_matrix.h5")


obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- CreateSeuratObject(counts = obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, project = "obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined", min.cells = 3, min.features = 200)
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined

dense.size <- object.size(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined)
dense.size

sparse.size <- object.size(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined)
sparse.size

obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined[["percent.mt"]] <- PercentageFeatureSet(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, pattern = "^MT-|^Mt-")

# Visualize QC metrics as a violin plot
jpeg(file="obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedPercentageFeatureSet.jpeg")
VlnPlot(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedFeatureScatter.jpeg")
plot1 <- FeatureScatter(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- NormalizeData(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, normalization.method = "LogNormalize", scale.factor = 10000)

obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- FindVariableFeatures(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, selection.method = "vst", nfeatures = 32000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined), 10)
write.csv(top10,"obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedtophighvariablegenes.csv")

top5 <- head(VariableFeatures(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined), 5)
write.csv(top5,"obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedtophighvariablegenes5.csv")

# plot variable features with and without labels
jpeg(file="obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedvariable_plot1.jpeg")
plot1 <- VariableFeaturePlot(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# plot variable features with and without labels
jpeg(file="obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedvariable5_plot1.jpeg")
plot1 <- VariableFeaturePlot(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined)
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- ScaleData(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, features = all.genes)

obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- RunPCA(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, features = VariableFeatures(object = obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined))

jpeg(file="obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedvizdim3load.jpeg")
VizDimLoadings(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, dims = 1:2, reduction = "pca")
dev.off()

obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- FindNeighbors(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, dims = 2:10)
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- FindClusters(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, resolution = 0.1)
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- RunUMAP(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, dims = 1:10)

jpeg(file="subset_lens_combined0.110dims.jpeg")
DimPlot(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, group.by	= "subset",reduction = "umap")
dev.off()

jpeg(file="umap_lens_combined0.110dims.jpeg")
DimPlot(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined,reduction = "umap")
dev.off()

obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- FindClusters(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, resolution = 0.1)
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- RunUMAP(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, dims = 2:10)
jpeg(file="obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined0.310dims.jpeg")
DimPlot(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, reduction = "umap")
dev.off()

obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- FindClusters(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, resolution = 0.4)
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined <- RunUMAP(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, dims = 1:10)
jpeg(file="obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined0.410dims.jpeg")
DimPlot(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, reduction = "umap", label=T)
dev.off()

# find all markers of cluster 2
cluster2.markers <- FindMarkers(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
write.csv(cluster2.markers, "obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combinedmarkers.csv")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined.markers <- FindAllMarkers(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(obj_integ_cornealimbus_combined_rs_0.1_ds_30_cornealimbus_combined.markers, "lens_combinedclustermarkers.csv")
