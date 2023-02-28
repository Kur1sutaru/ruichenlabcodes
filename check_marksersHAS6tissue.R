# library
source("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/R_functions_scATAC.R")
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(reticulate)

#user parameters
objectA   <- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scRNA_HAS/integration/predict_Robjs_rna/out_integration_ciliary_combined/obj_integ_ciliary_combined_rs_0.1_ds_30_ciliary_combined.rds")
markerSet <- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/marker_genes_HAS/test_checkMarkers/obj_marker_HAS_list_6T.rds")
prefix <- "cilia"
dirOut <- getwd()

# check markers for 6 tissue types  
checkMarkersHsHas(objectA,markerSet,prefix,dirOut)


