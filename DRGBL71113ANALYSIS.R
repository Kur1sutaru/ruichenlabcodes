## quality control sample bl-71113
setwd("")
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(reticulate)
library(dplyr)
library(patchwork)
library(ggplot2)
library(dplyr)

bl_71113 <- Read10X_h5(filename = "D:/Baylor/DRG/bl-71113/filtered_feature_bc_matrix.h5")
bl_71113 <- CreateSeuratObject(counts = bl_71113, project = "bl_71113", min.cells = 3, min.features = 200)
bl_71113
bl_71113[["percent.mt"]] <- PercentageFeatureSet(bl_71113, pattern = "^MT-|^Mt-")

# 5% for single nuclei, 10% for single cell
bl_71113 <- subset(bl_71113, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Visualize QC metrics as a violin plot
jpeg(file="bl_71113PercentageFeatureSet.jpeg")
VlnPlot(bl_71113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

jpeg(file="bl_71113FeatureScatter.jpeg")
plot1 <- FeatureScatter(bl_71113, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bl_71113, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()


bl_71113 <- NormalizeData(bl_71113, normalization.method = "LogNormalize", scale.factor = 10000)

bl_71113 <- FindVariableFeatures(bl_71113, selection.method = "vst", nfeatures = 21525)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bl_71113), 10)
write.csv(top10,"bl_71113tophighvariablegenes.csv")

top5 <- head(VariableFeatures(bl_71113), 5)
write.csv(top5,"bl_71113tophighvariablegenes5.csv")

# plot variable features with and without labels
jpeg(file="bl_71113variable_plot1.jpeg")
plot1 <- VariableFeaturePlot(bl_71113)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# plot variable features with and without labels
jpeg(file="bl_71113variable5_plot1.jpeg")
plot1 <- VariableFeaturePlot(bl_71113)
plot2 <- LabelPoints(plot = plot1, points = top5, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(bl_71113)
bl_71113 <- ScaleData(bl_71113, features = all.genes)

bl_71113 <- RunPCA(bl_71113, features = VariableFeatures(object = bl_71113))

jpeg(file="bl_71113vizdim3load.jpeg")
VizDimLoadings(bl_71113, dims = 1:2, reduction = "pca")
dev.off()

bl_71113 <- FindNeighbors(bl_71113, dims = 1:30)
bl_71113 <- FindClusters(bl_71113, resolution = 0.8)
bl_71113 <- RunUMAP(bl_71113, dims = 1:30)

DimPlot(bl_71113,reduction = "umap")

VlnPlot(bl_71113, features = c("Atf3", "Calca", "Trpm8", "Scn7a", "Dcn"))


jpeg(file="umap_sample66091res008.jpeg")
DimPlot(bl_71113,reduction = "umap")
dev.off()

saveRDS(bl_71113, "bl_71113.rds")
################################ scpred prediction #####################################
library("scPred")
library("glue")
library("dplyr")
library("Seurat")
library("magrittr")
library("caret")
options(Seurat.object.assay.version = "v5")
## to avoid compatibility errors, install by devtools::install_github(repo="powellgenomicslab/scPred",  ref="9f407b7436f40d44224a5976a94cc6815c6e837f")
#Training the model
DRG_neurons <- DRG_neurons %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

## PLAN B
all.genes <- rownames(DRG_neurons)
DRG_neurons <- ScaleData(DRG_neurons, features = all.genes)
DRG_neurons <- RunPCA(DRG_neurons, features = VariableFeatures(object = DRG_neurons))
DRG_neurons <- RunUMAP(DRG_neurons, dims = 1:30)
DimPlot(DRG_neurons, group.by = "Atlas_annotation", label = TRUE, repel = TRUE)
DRG_neurons <- getFeatureSpace(DRG_neurons, "Atlas_annotation")
DRG_neurons <- trainModel(DRG_neurons)
get_probabilities(DRG_neurons) %>% head()
get_scpred(DRG_neurons)

#plot_probabilities(DRG_neurons)


### After training the model with the DRG_neurons, lets predict the query - our data set
# Normalize the data
#bl_71113   <- NormalizeData(bl_71113)

# Find variable features
bl_71113   <- FindVariableFeatures(bl_71113)

# Scale the data
bl_71113   <- ScaleData(bl_71113)

# Perform PCA
bl_71113   <- RunPCA(bl_71113)

# Check available assays
Assays(bl_71113)

# Make sure "integrated" is the default assay for predictions
bl_71113 <- scPredict(bl_71113, DRG_neurons)
DimPlot(bl_71113, group.by = "scpred_prediction", reduction = "scpred")
bl_71113   <- RunUMAP(bl_71113, reduction = "scpred", dims = 1:30)
DimPlot(bl_71113, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
# Store the trained scPred model in the 'misc' slot of the DRG_neurons Seurat object
bl_71113@misc$scPred_model <- DRG_neurons@tools$scPred
# Calculate proportions of cell types
cell_proportions <- table(bl_71113@meta.data$scpred_prediction)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)
# Load ggplot2 for plotting
library(ggplot2)

# Assuming you have already created the data frame with cell type proportions:
cell_proportions <- table(bl_71113@meta.data$scpred_prediction)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)

# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add the cell counts above the bars
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from scPred drg_bl71113") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


saveRDS(bl_71113, "bl_71113_annotatedneurons.rds")



### Accessing classifiers
crossTab(bl_71113, "Atlas_annotation", "scpred_prediction")
crossTab(bl_71113, "Atlas_annotation", "scpred_prediction", output = "prop")
get_classifiers(DRG_neurons)

# Each model can be normally treated using the caret enviroment. 
# For example, we can plot the performance resamples using the plot.train:
caret::plot.train(get_classifiers(DRG_neurons)[["Chondrocytes"]])



### non neurons
## PLAN B
all.genes <- rownames(DRG_nonneurons)
DRG_nonneurons <- ScaleData(DRG_nonneurons, features = all.genes)
DRG_nonneurons <- RunPCA(DRG_nonneurons, features = VariableFeatures(object = DRG_nonneurons))
DRG_nonneurons <- RunUMAP(DRG_nonneurons, dims = 1:30)
DimPlot(DRG_nonneurons, group.by = "Atlas_annotation", label = TRUE, repel = TRUE)
DRG_nonneurons <- getFeatureSpace(DRG_nonneurons, "Atlas_annotation")
DRG_nonneurons <- trainModel(DRG_nonneurons)
get_probabilities(DRG_nonneurons) %>% head()
get_scpred(DRG_nonneurons)
bl_71113_annotated <- scPredict(bl_71113_annotated, DRG_nonneurons)
DimPlot(bl_71113_annotated, group.by = "scpred_prediction", reduction = "scpred")
bl_71113_annotated   <- RunUMAP(bl_71113_annotated, reduction = "scpred", dims = 1:30)
DimPlot(bl_71113_annotated, group.by = "scpred_prediction", label = TRUE, repel = TRUE)

# Store the trained scPred model in the 'misc' slot of the DRG_neurons Seurat object
bl_71113_annotated@misc$scPred_model <- DRG_nonneurons@tools$scPred

# Assuming you have already created the data frame with cell type proportions:
cell_proportions <- table(bl_71113_annotated@meta.data$scpred_prediction)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)


# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add the cell counts above the bars
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from scPred drg_bl71113 non neuronal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


### Subset unassigned cells and re cluster them
# View unique predictions in the scpred_prediction column
unique(bl_71113_annotated@meta.data[["scpred_prediction"]])
# Subset the cells where scpred_prediction is "unassigned"
DRG_nonneuronal <- subset(bl_71113_annotated, subset = scpred_prediction == "unassigned")
# Check the number of cells in the unassigned cluster
table(DRG_nonneuronal@meta.data[["scpred_prediction"]])
DRG_nonneuronal[["percent.mt"]] <- PercentageFeatureSet(DRG_nonneuronal, pattern = "^MT-|^Mt-")

# 5% for single nuclei, 10% for single cell
DRG_nonneuronal <- subset(DRG_nonneuronal, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Visualize QC metrics as a violin plot
jpeg(file="DRG_nonneuronalPercentageFeatureSet.jpeg")
VlnPlot(DRG_nonneuronal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
DRG_nonneuronal <- NormalizeData(DRG_nonneuronal, normalization.method = "LogNormalize", scale.factor = 10000)

DRG_nonneuronal <- FindVariableFeatures(DRG_nonneuronal, selection.method = "vst", nfeatures = 21525)
all.genes <- rownames(DRG_nonneuronal)
DRG_nonneuronal <- ScaleData(DRG_nonneuronal, features = all.genes)

DRG_nonneuronal <- RunPCA(DRG_nonneuronal, features = VariableFeatures(object = DRG_nonneuronal))

jpeg(file="DRG_nonneuronalvizdim3load.jpeg")
VizDimLoadings(DRG_nonneuronal, dims = 1:2, reduction = "pca")
dev.off()

DRG_nonneuronal <- FindNeighbors(DRG_nonneuronal, dims = 1:30)
DRG_nonneuronal <- FindClusters(DRG_nonneuronal, resolution = 0.8)
DRG_nonneuronal <- RunUMAP(DRG_nonneuronal, dims = 1:30)

DimPlot(DRG_nonneuronal,reduction = "umap")



### gene expression dot plot to annotate non neuronal clusters
# List of genes to analyze
  genes_of_interest <- c("Ednrb", "Fabp7", "Cadm2", "Scn7a", "Cldn5", "Ptgds","Mgp", "Dcn", "Egfl7", "Pecam1", "Prx", "Lyz2", "Mrc1","Mbp", "S100a8", "Olig1", "Kcnj8", "Mpz", "Cd3e")

# Generate the dot plot
  DotPlot(DRG_nonneuronal, features = genes_of_interest) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Gene Expression DotPlot", x = "Gene", y = "Cluster")

#violin plot Satglia
VlnPlot(DRG_nonneuronal, features = c("Ednrb", "Fabp7"))
#violin plot Schwann_N
VlnPlot(DRG_nonneuronal, features = c("Cadm2", "Scn7a"))
#violin plot Fibroblast
VlnPlot(DRG_nonneuronal, features = c("Ptgds", "Mgp", "Dcn"))
#violin plot Endothelial
VlnPlot(DRG_nonneuronal, features = c("Egfl7", "Pecam1", "Cldn5"))
#violin plot Schwann_M
VlnPlot(DRG_nonneuronal, features = c("Prx", "Mbp", "Mpz"))
#violin plot Immune
VlnPlot(DRG_nonneuronal, features = c("Mrc1", "Cd4", "Lyz2"))
#violin plot Pericyte
VlnPlot(DRG_nonneuronal, features = c("Des", "Myh11", "Kcnj8", "Notch3"))


DimPlot(DRG_neurons, reduction = "umap", group.by = "Atlas_annotation", label = TRUE, pt.size = 0.5)
# Assuming you have already created the data frame with cell type proportions:
cell_proportions <- table(DRG_neurons@meta.data$Atlas_annotation)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)


# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add the cell counts above the bars
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from DRG_neurons meta atlas") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


### new clusters ids
new.cluster.ids <- c("Fibroblast1", "Pericyte", "Satglia", "Schwann_N", "Schwann_M", "Schwann_M",
                     "Endothelial", "Schwann_N", "Schwann_N", "Schwann_M", "Fibroblast", "Immune", "Fibroblast3")
names(new.cluster.ids) <- levels(DRG_nonneuronal)
DRG_nonneuronal <- RenameIdents(DRG_nonneuronal, new.cluster.ids)
DimPlot(DRG_nonneuronal, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Assuming you have already created the data frame with cell type proportions:
cell_proportions <- table(DRG_nonneuronal@meta.data$Atlas_annotation)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)


# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add the cell counts above the bars
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from DRG_neurons meta atlas") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# View available metadata columns
colnames(DRG_nonneuronal@meta.data)

# Remove multiple columns
DRG_nonneuronal@meta.data <- DRG_nonneuronal@meta.data[, !(colnames(bl_71113_annotated@meta.data) %in% c("scpred_Mrgprd", "scpred_Trpm8", "scpred_Rxfp1", "scpred_Mrgpra3_plusTrpv1", "scpred_Pvalb", "scpred_Ntrk3high_plusS100a16", "scpred_Calca_plusAdra2a", "scpred_Calca_plusSstr2", "scpred_Calca_plusBmpr1b", "scpred_Mrgpra3_plusMrgprb4", "scpred_Ntrk3low_plusNtrk2", "scpred_Atf3", "scpred_no_rejection", "scpred_Th", "scpred_Calca_plusDcn", "scpred_Calca_plusSmr2", "scpred_Sst", "scpred_Calca_plusOprk1", "scpred_Ntrk3high_plusNtrk2", "scpred_max"))]

# Verify the columns have been removed
colnames(DRG_nonneuronal@meta.data)
# Copy the active identity to a new metadata column named "Atlas_annotation"
DRG_nonneuronal@meta.data$Atlas_annotation <- Idents(DRG_nonneuronal)

# Verify that the column has been added
head(DRG_nonneuronal@meta.data$Atlas_annotation)
# Assuming you have already created the data frame with cell type proportions:
cell_proportions <- table(DRG_nonneuronal@meta.data$Atlas_annotation)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)


# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add the cell counts above the bars
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from DRG non neuronal sample bl71113") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

## cell proportion non neuronal meta atlas
# Assuming you have already created the data frame with cell type proportions:
cell_proportions <- table(DRG_nonneurons@meta.data$Atlas_annotation)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)


# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +  # Add the cell counts above the bars
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from DRG non neuronal meta atlas") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
saveRDS(DRG_nonneuronal, "DRG_nonneuronalbl71113.rds")
