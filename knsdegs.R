setwd("D:/baylor/KNS")
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)
library(mclust)
# load data

##Quality control

K_B6_NSDS5[["percent.mt"]] <- PercentageFeatureSet(K_B6_NSDS5, pattern = "^MT-|^Mt-")

# Visualize QC metrics as a violin plot
jpeg(file="K_B6_NSDS5PercentageFeatureSet.jpeg")
VlnPlot(K_B6_NSDS5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="K_B6_NSDS5FeatureScatter.jpeg")
plot1 <- FeatureScatter(K_B6_NSDS5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(K_B6_NSDS5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()



###Normalization and pre clustering
K_B6_NSDS5 <- NormalizeData(K_B6_NSDS5, normalization.method = "LogNormalize", scale.factor = 10000)

K_B6_NSDS5 <- FindVariableFeatures(K_B6_NSDS5, selection.method = "vst", nfeatures = 36601)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(K_B6_NSDS5), 10)
write.csv(top10,"K_B6_NSDS5tophighvariablegenes.csv")

top5 <- head(VariableFeatures(K_B6_NSDS5), 5)
write.csv(top5,"K_B6_NSDS5tophighvariablegenes5.csv")

# plot variable features with and without labels
jpeg(file="K_B6_NSDS5variable_plot1.jpeg")
plot1 <- VariableFeaturePlot(K_B6_NSDS5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# plot variable features with and without labels
jpeg(file="K_B6_NSDS5variable5_plot1.jpeg")
plot1 <- VariableFeaturePlot(K_B6_NSDS5)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(K_B6_NSDS5)
K_B6_NSDS5 <- ScaleData(K_B6_NSDS5, features = all.genes)

K_B6_NSDS5 <- RunPCA(K_B6_NSDS5, features = VariableFeatures(object = K_B6_NSDS5))

jpeg(file="K_B6_NSDS5vizdim3load.jpeg")
VizDimLoadings(K_B6_NSDS5, dims = 1:2, reduction = "pca")
dev.off()

K_B6_NSDS5 <- FindNeighbors(K_B6_NSDS5, dims = 1:30)
K_B6_NSDS5 <- FindClusters(K_B6_NSDS5, resolution = 0.5)
K_B6_NSDS5 <- RunUMAP(K_B6_NSDS5, dims = 1:30)

#To generate umaps for each sample in orig.ident
jpeg(file="splitcellltypesSTROMALorigidentK_B6_NSDS5vizdim3load.jpeg")
DimPlot(K_B6_NSDS5_Mins_Stromal, group.by = "celltype", split.by = "orig.ident", reduction = "umap")
dev.off()




# find markers for every cluster compared to all remaining cells, report only the positive - generate the Differentially expressed genes - DEGs list
# ones
K_B6_NSDS5_Mins_Stromal <- FindAllMarkers(K_B6_NSDS5_Mins_Stromal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
K_B6_NSDS5_Mins_Stromal %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(K_B6_NSDS5_Mins_Stromal, "degsK_B6_NSDS5_Mins_Stromal.csv")

SaveH5Seurat(K_B6_NSDS5, overwrite = TRUE)
SaveH5Seurat(K_B6_NSDS5, filename = "K_B6_NSDS5.h5Seurat")
Convert("K_B6_NSDS5.h5Seurat", dest = "h5ad")
saveRDS(K_B6_NSDS5, "K_B6_NSDS5.rds")

# Using group.by to make the comparison between NS and DS5
testdegsK_B6_NSDS5 <- FindMarkers(
  K_B6_NSDS5_Mins_Stromal,
  group_by = "celltypes",
  split.by = "orig.ident",
  ident.1 = "PF_61353",
  ident.2 = "PF_63067")


markersboth <- FindMarkers(K_B6_NSDS5_Mins_Stromal, ident.1 = "PF_61353", ident.2 = "PF_63067", group.by = 'orig.ident')
K_B6_NSDS5_Mins_Stromal %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(K_B6_NSDS5_Mins_Stromal, "splitdegsK_B6_NSDS5_Mins_Stromal.csv")


#Aditional steps - optional
# Load the Seurat package and your Seurat object
library(Seurat)
my_seurat <- Read10X(data.dir = "path/to/data")
my_seurat <- CreateSeuratObject(counts = my_seurat$`Gene Expression`)

# Identify cell types using clustering and/or manual annotation
# For example:
my_seurat <- FindNeighbors(my_seurat, dims = 1:10)
my_seurat <- FindClusters(my_seurat, resolution = 0.5)
K_B6_NSDS5$celltype <- Idents(K_B6_NSDS5)

# Create the cell type proportion barplot
VlnPlot(K_B6_NSDS5, features = "celltype", pt.size = 0, group.by = "celltype")
VlnPlot(K_B6_NSDS5, features = "celltype", pt.size = 0, group.by = "celltype", geom = "bar")

#Feature expression heatmap

jpeg(file="objB6TGNuheatmap.jpeg")
DoHeatmap(
  K_B6_NSDS5,
  features = NULL,
  cells = NULL,
  group.by = "celltype",
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

counts <- GetAssyData(K_B6_NSDS5, assay="RNA", slot="data")
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


#Cell proportion barplot
pt <- table(Idents(K_B6_NSDS5), K_B6_NSDS5$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("ON integrated") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

#Other cell proportion barplot

# Extract the cell type or cluster information from the Seurat object
# Replace 'ClusterName' with the actual metadata column containing cell type or cluster labels
cell_type_info <- onsamplesoctober$ClusterName

# Calculate the proportions
cell_type_counts <- table(cell_type_info)
proportions <- prop.table(cell_type_counts) * 100  # Calculate proportions and multiply by 100 for percentages

# Convert proportions to a data frame for ggplot2
cell_type_df <- data.frame(CellType = names(proportions), Proportion = proportions)

# Create the barplot with ggplot2
cell_proportion_plot <- ggplot(cell_type_df, aes(x = CellType, y = Proportion)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Cell Type Proportions",
       x = "Cell Type",
       y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

# Display the plot
print(cell_proportion_plot)


#violin and feature plots

VlnPlot(K_B6_NSDS5, features = c("BEST1","RPE65"))
FeaturePlot(K_B6_NSDS5, features = c("SOX9","OTX2"))


saveRDS(K_B6_NSDS5, file = "K_B6_NSDS5.rds")


