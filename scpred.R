setwd()
library("scPred")
library("glue")
library("dplyr")
library("Seurat")
library("magrittr")
library("caret")
options(Seurat.object.assay.version = "v5")
## to avoid compatibility errors, install by devtools::install_github(repo="powellgenomicslab/scPred",  ref="9f407b7436f40d44224a5976a94cc6815c6e837f")
#Training the model
reference <- scPred::longbone_stromal 
query <- scPred::pbmc_2
reference <- longbone_stromal %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE)
reference <- getFeatureSpace(reference, "cell_type")
reference <- trainModel(reference)
get_probabilities(reference) %>% head()
get_scpred(reference)
plot_probabilities(reference)

### After training the model with the reference, lets predict the query - our data set
# Normalize the data
kolb66092subset   <- NormalizeData(kolb66092subset  )

# Find variable features
kolb66092subset   <- FindVariableFeatures(kolb66092subset  )

# Scale the data
kolb66092subset   <- ScaleData(kolb66092subset  )

# Perform PCA
kolb66092subset   <- RunPCA(kolb66092subset  )

# Check available assays
Assays(kolb66092subset  )

# Make sure "integrated" is the default assay for predictions
kolb66092subset   <- scPredict(kolb66092subset  , reference)
DimPlot(kolb66092subset  , group.by = "scpred_prediction", reduction = "scpred")
kolb66092subset   <- RunUMAP(kolb66092subset  , reduction = "scpred", dims = 1:30)
DimPlot(kolb66092subset  , group.by = "scpred_prediction", label = TRUE, repel = TRUE)
# Store the trained scPred model in the 'misc' slot of the reference Seurat object
kolb66092subset  @misc$scPred_model <- reference@tools$scPred
# Calculate proportions of cell types
cell_proportions <- table(kolb66092subset  @meta.data$scpred_prediction)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)
# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from scPred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

saveRDS(wt66089subset, "wtlb66089subset_annotated.rds")
saveRDS(kolb66092subset , "kolb66092subset _annotated.rds")
saveRDS(kolb66092subset , "kolb66092subset _annotated.rds")
saveRDS(kolb66092subset, "kolb66092subset_annotated.rds")


### Accessing classifiers
crossTab(kolb66092subset  , "cell_type", "scpred_prediction")
crossTab(kolb66092subset  , "cell_type", "scpred_prediction", output = "prop")
get_classifiers(reference)

# Each model can be normally treated using the caret enviroment. 
# For example, we can plot the performance resamples using the plot.train:
caret::plot.train(get_classifiers(reference)[["Chondrocytes"]])


# Calculate proportions of cell types
cell_proportions <- table(reference@meta.data$cell_type)
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("CellType", "Count")
cell_proportions_df$Proportion <- cell_proportions_df$Count / sum(cell_proportions_df$Count)
# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot
ggplot(cell_proportions_df, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Predicted Cell Type", y = "Proportion", title = "Cell Type Proportions from scPred") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



