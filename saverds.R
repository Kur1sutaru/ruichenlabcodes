library(Seurat)
library(SeuratDisk)
library(SeuratObject)

## To save h5 or h5ad files into .rds
b6tgnu<- Read10X_h5(filename = "D:/baylor/TRIGEMINAL/B6_TG_Nu/b6tgnu.h5")
# Save an object to a file
saveRDS(b6tgnu, file = "b6tgnu.rds")
