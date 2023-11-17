library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Seurat)
library(ggplot2)


#### https://satijalab.org/seurat/archive/v3.1/immune_alignment.html

DimPlot(K_B6_NSDS5_Mins_Stromal, reduction = "umap", split.by = "orig.ident")

#To generate degs each clusters - ident = cluster number

DefaultAssay(K_B6_NSDS5_Mins_Stromal) <- "RNA"
cluster1.markers <- FindConservedMarkers(K_B6_NSDS5_Mins_Stromal, ident.1 = "Mast", grouping.var = "orig.ident", verbose = FALSE)
write.csv(cluster1.markers, "Mastdegs.csv")
head(cluster1.markers)

