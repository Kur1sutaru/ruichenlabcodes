# Load required library
library(Seurat)

# Assuming 'seurat_obj' is your Seurat object

# Extract counts from the 'RNA' assay
bulk_counts <- GetAssayData(object = seurat_obj, slot = "counts")

# If you want to get counts for specific genes, you can subset the count matrix
# For example, to extract counts for genes "Gene1", "Gene2", and "Gene3"
selected_genes <- c("Gene1", "Gene2", "Gene3")
bulk_counts_selected <- bulk_counts[selected_genes, ]

# Now 'bulk_counts' contains the count matrix
# Load required library
library(Seurat)

# Assuming 'seurat_obj' is your Seurat object

# Check available assays in the Seurat object
Assays(seurat_obj)

# If "RNA" is not already the active assay, set it as the active assay
seurat_obj <- SetAssay(seurat_obj, assay = "RNA")

# Now the "RNA" assay is set as the active assay in the Seurat object
