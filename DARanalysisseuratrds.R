library(Seurat)

# Load the Seurat object containing the ATAC-seq data
atac_seq_data <- readRDS("atac_seq_data.rds")

# Normalize the ATAC-seq signal across cells
atac_seq_data <- NormalizeData(atac_seq_data, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the highly variable regions
atac_seq_data <- FindVariableFeatures(atac_seq_data, selection.method = "vst", nfeatures = 2000)

# Scale the ATAC-seq signal
atac_seq_data <- ScaleData(atac_seq_data)

# Get the conditions
conditions <- unique(atac_seq_data@meta.data$condition)

# Perform differential analysis to identify differentially accessible regions (DARs)
dars <- FindMarkers(atac_seq_data, ident.1 = conditions[1], ident.2 = conditions[2], min.pct = 0.25, logfc.threshold = 0.25)

# Store the differentially accessible regions (DARs)
saveRDS(dars, "dars.rds")
