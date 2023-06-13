setwd("D:/baylor/HumanAnteriorSegment/fetalHAS")
library(Seurat)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)
library(mclust)
# load data

objB6TGNu <- Read10X_h5(filename = "D:/baylor/TRIGEMINAL/objB6TGNu/analysis/input/filtered_feature_bc_matrix.h5")


objB6TGNu <- CreateSeuratObject(counts = objB6TGNu, project = "objB6TGNu", min.cells = 3, min.features = 200)
objB6TGNu

dense.size <- object.size(objB6TGNu)
dense.size

sparse.size <- object.size(objB6TGNu)
sparse.size

objB6TGNu[["percent.mt"]] <- PercentageFeatureSet(objB6TGNu, pattern = "^MT-|^Mt-")

# Visualize QC metrics as a violin plot
jpeg(file="objB6TGNuPercentageFeatureSet.jpeg")
VlnPlot(objB6TGNu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="objB6TGNuFeatureScatter.jpeg")
plot1 <- FeatureScatter(objB6TGNu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(objB6TGNu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


objB6TGNu <- NormalizeData(objB6TGNu, normalization.method = "LogNormalize", scale.factor = 10000)

objB6TGNu <- FindVariableFeatures(objB6TGNu, selection.method = "vst", nfeatures = 32000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(objB6TGNu), 10)
write.csv(top10,"objB6TGNutophighvariablegenes.csv")

top5 <- head(VariableFeatures(objB6TGNu), 5)
write.csv(top5,"objB6TGNutophighvariablegenes5.csv")

# plot variable features with and without labels
jpeg(file="objB6TGNuvariable_plot1.jpeg")
plot1 <- VariableFeaturePlot(objB6TGNu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# plot variable features with and without labels
jpeg(file="objB6TGNuvariable5_plot1.jpeg")
plot1 <- VariableFeaturePlot(objB6TGNu)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(objB6TGNu)
objB6TGNu <- ScaleData(objB6TGNu, features = all.genes)

objB6TGNu <- RunPCA(objB6TGNu, features = VariableFeatures(object = objB6TGNu))

jpeg(file="objB6TGNuvizdim3load.jpeg")
VizDimLoadings(objB6TGNu, dims = 1:2, reduction = "pca")
dev.off()

objB6TGNu <- FindNeighbors(objB6TGNu, dims = 1:10)
objB6TGNu <- FindClusters(objB6TGNu, resolution = 0.1)
objB6TGNu <- RunUMAP(objB6TGNu, dims = 1:10)

jpeg(file="subset_lens_combined0.110dims.jpeg")
DimPlot(objB6TGNu, group.by	= "subset",reduction = "umap")
dev.off()

#Feature expression heatmap

jpeg(file="objB6TGNuheatmap.jpeg")
DoHeatmap(
  objB6TGNu,
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

#Heatmap with a given gene list

counts <- GetAssyData(objB6TGNu, assay="RNA", slot="data")
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