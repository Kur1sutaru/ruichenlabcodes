setwd("D:/baylor/DRG atlas")
library(Seurat)
library(ggplot2)

# To verify how many cell types we have to do DEG analysis
levels(drg.combined)

# [1] "Putative C-LTMRs"             "Pruritogen receptor enriched" "Putative silent nociceptors"  "TRPA1+ nociceptors"          
# [5] "PENK+ nociceptors"            "Ad HTMRs"                     "Cold nociceptors"             "Aß nociceptors"              
# [9] "Ad LTMRs"                     "Aß RA LTMRs"                  "Aß SA LTMRs"                  "Proprioceptors"

#Find DEGs between Cold nociceptors and Aß nociceptors
nociceptormarkers <- FindMarkers(drg.combined, ident.1 = "Cold nociceptors", ident.2 = "Aß nociceptors")
# view results
head(nociceptormarkers)

# Save results
write.csv(nociceptormarkers, "coldnociceptorsxabnociceptors_DEG.csv")

#Find DEGs between Putative C-LTMRs and PENK+ nociceptors
nociceptormarkers2 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "PENK+ nociceptors")
write.csv(nociceptormarkers2, "C-LTMRsxPENKnoci_DEG.csv")

#Find DEGs between Putative C-LTMRs and Ad LTMRs
nociceptormarkers3 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Ad LTMRs")
write.csv(nociceptormarkers3, "C-LTMRsxAdLTMRs_DEG.csv")

#Find DEGs between Putative C-LTMRs and Pruritogen receptor enriched 
nociceptormarkers4 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Pruritogen receptor enriched")
write.csv(nociceptormarkers4, "C-LTMRsxpluritogenreceptors_DEG.csv")

#Find DEGs between Putative C-LTMRs and Ad HTMRs
nociceptormarkers5 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Ad HTMRs")
write.csv(nociceptormarkers5, "C-LTMRsxAdHTMRs_DEG.csv")

#Find DEGs between Putative C-LTMRs and Aß RA LTMRs
nociceptormarkers5 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Aß RA LTMRs")
write.csv(nociceptormarkers5, "C-LTMRsxAßRALTMRs_DEG.csv")

#Find DEGs between Putative C-LTMRs and Putative silent nociceptors
nociceptormarkers6 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Putative silent nociceptors")
write.csv(nociceptormarkers6, "C-LTMRsxPutativesilnocicep_DEG.csv")

#Find DEGs between Putative C-LTMRs and Putative silent nociceptors
nociceptormarkers7 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Cold nociceptors")
write.csv(nociceptormarkers7, "C-LTMRsxcoldnocicep_DEG.csv")

#Find DEGs between Putative C-LTMRs and Putative silent nociceptors
nociceptormarkers8 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Aß SA LTMRs")
write.csv(nociceptormarkers8, "C-LTMRsxAßSALTMRs_DEG.csv")

#Find DEGs between Putative C-LTMRs and Putative silent nociceptors
nociceptormarkers9 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "TRPA1+ nociceptors")
write.csv(nociceptormarkers9, "C-LTMRsxTRPA1nocis_DEG.csv")

#Find DEGs between Putative C-LTMRs and Putative silent nociceptors
nociceptormarkers10 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Aß nociceptors")
write.csv(nociceptormarkers10, "C-LTMRsxAßnocis_DEG.csv")

#Find DEGs between Putative C-LTMRs and Putative silent nociceptors
nociceptormarkers11 <- FindMarkers(drg.combined, ident.1 = "Putative C-LTMRs", ident.2 = "Proprioceptors")
write.csv(nociceptormarkers11, "C-LTMRsxProprioceptors_DEG.csv")

# Single cell heatmap of feature expression
DoHeatmap(subset(drg.combined, downsample = 100), features = drg.combined@active.ident, size = 3)

# Plot UMAP, coloring cells by cell type (currently stored in object@ident)
DimPlot(drg.combined, reduction = "umap")

data("drg.combined")
DoHeatmap(object = drg.combined)





