setwd("D:/baylor/TRIGEMINAL/B6_TG_Nu")
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)
library(mclust)
# load data

B6_TG_Nu <- Read10X_h5(filename = "D:/baylor/TRIGEMINAL/B6_TG_Nu/b6tgnu.h5")
B6_TG_Nu <-LoadH5Seurat("D:/baylor/TRIGEMINAL/B6_TG_Nu/Neuron_paper_cell_label/objintegrated.h5seurat")

B6_TG_Nu <- CreateSeuratObject(counts = B6_TG_Nu, project = "B6_TG_Nu", min.cells = 3, min.features = 200)
B6_TG_Nu

dense.size <- object.size(B6_TG_Nu)
dense.size

sparse.size <- object.size(B6_TG_Nu)
sparse.size

B6_TG_Nu[["percent.mt"]] <- PercentageFeatureSet(B6_TG_Nu, pattern = "^MT-|^Mt-")

# Visualize QC metrics as a violin plot
jpeg(file="B6_TG_NuPercentageFeatureSet.jpeg")
VlnPlot(B6_TG_Nu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="B6_TG_NuFeatureScatter.jpeg")
plot1 <- FeatureScatter(B6_TG_Nu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(B6_TG_Nu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


B6_TG_Nu <- NormalizeData(B6_TG_Nu, normalization.method = "LogNormalize", scale.factor = 10000)

B6_TG_Nu <- FindVariableFeatures(B6_TG_Nu, selection.method = "vst", nfeatures = 32000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(B6_TG_Nu), 10)
write.csv(top10,"B6_TG_Nutophighvariablegenes.csv")

top5 <- head(VariableFeatures(B6_TG_Nu), 5)
write.csv(top5,"B6_TG_Nutophighvariablegenes5.csv")

# plot variable features with and without labels
jpeg(file="B6_TG_Nuvariable_plot1.jpeg")
plot1 <- VariableFeaturePlot(B6_TG_Nu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# plot variable features with and without labels
jpeg(file="B6_TG_Nuvariable5_plot1.jpeg")
plot1 <- VariableFeaturePlot(B6_TG_Nu)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(B6_TG_Nu)
B6_TG_Nu <- ScaleData(B6_TG_Nu, features = all.genes)

B6_TG_Nu <- RunPCA(B6_TG_Nu, features = VariableFeatures(object = B6_TG_Nu))

jpeg(file="B6_TG_Nuvizdim3load.jpeg")
VizDimLoadings(B6_TG_Nu, dims = 1:2, reduction = "pca")
dev.off()

B6_TG_Nu <- FindNeighbors(B6_TG_Nu, dims = 1:10)
B6_TG_Nu <- FindClusters(B6_TG_Nu, resolution = 0.1)
B6_TG_Nu <- RunUMAP(B6_TG_Nu, dims = 1:10)

jpeg(file="subset_lens_combined0.110dims.jpeg")
DimPlot(B6_TG_Nu, group.by = "name",reduction = "umap")
dev.off()

jpeg(file="umap_lens_combined0.110dims.jpeg")
DimPlot(B6_TG_Nu,reduction = "umap")
dev.off()

B6_TG_Nu <- FindClusters(B6_TG_Nu, resolution = 0.1)
B6_TG_Nu <- RunUMAP(B6_TG_Nu, dims = 1:10)
jpeg(file="B6_TG_Nu0.310dims.jpeg")
DimPlot(B6_TG_Nu, reduction = "umap")
dev.off()

B6_TG_Nu <- FindClusters(B6_TG_Nu, resolution = 0.4)
B6_TG_Nu <- RunUMAP(B6_TG_Nu, dims = 1:10)
jpeg(file="B6_TG_Nu0.410dims.jpeg")
DimPlot(B6_TG_Nu, reduction = "umap", label=T)
dev.off()

# find all markers of cluster 2
cluster2.markers <- FindMarkers(B6_TG_Nu, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
write.csv(cluster2.markers, "B6_TG_Numarkers.csv")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
B6_TG_Nu.markers <- FindAllMarkers(B6_TG_Nu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
B6_TG_Nu.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(B6_TG_Nu.markers, "tg_combinedclustermarkers.csv")



#Aditional steps
# Load the Seurat package and your Seurat object
library(Seurat)
my_seurat <- Read10X(data.dir = "path/to/data")
my_seurat <- CreateSeuratObject(counts = my_seurat$`Gene Expression`)

# Identify cell types using clustering and/or manual annotation
# For example:
my_seurat <- FindNeighbors(my_seurat, dims = 1:10)
my_seurat <- FindClusters(my_seurat, resolution = 0.5)
my_seurat$cell_type <- Idents(my_seurat)

# Create the cell type proportion barplot
VlnPlot(my_seurat, features = "cell_type", pt.size = 0, group.by = "cell_type")
VlnPlot(my_seurat, features = "cell_type", pt.size = 0, group.by = "cell_type", geom = "bar")

#Feature expression heatmap

jpeg(file="objB6TGNuheatmap.jpeg")
DoHeatmap(
  B6_TG_Nu,
  features = NULL,
  cells = NULL,
  group.by = "name",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
dev.off()

#Heatmap with a given gene list

counts <- GetAssyData(B6_TG_Nu, assay="RNA", slot="data")
genes <- c("Iba1","Gfap","Iba1","Atf3","cFos","Fos","Il1b","Il6","Nos2","Nox","Ccl2","Cd68","Itgam","Tac1","Calca","Trpm8","Trpv1","S100b","Gfra2",
           "Pou4F2","Gal","Cd55","Scn11a","Fxyd7","Ngfr","Nefh","Hapln4","Cbln2","Kcnab1","Sst","Il31ra","Apoe","Fabp7","Mpz","Gldn","Scn7a","Dcn","Pdgfra","Mgp","Alpl",
           "Cd74","Igfbp7","Tinagl1","Htr1f","Klf6","Klf9","Mt1","Egr1","Egr2","Cyr61","Nr4a1","Ctgf","Jag1","Lrp1")
counts <- as.matrix(counts[rownames(counts) %in% genes, ])

jpeg(file="genelistB6TGNuheatmap.jpeg")
DoHeatmap(
  counts,
  features = NULL,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)


dev.off()



saveRDS(B6_TG_Nu, file = "B6_TG_Nu.rds")
