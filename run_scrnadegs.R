### set dir -----
#source("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scRNA_HAS/scRNA_DEGs/run_scrnadegs.R")
source("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scRNA_HAS/scRNA_DEGs/run_scrnadegs.R")
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(reticulate)

objectA <- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scRNA_HAS/integration/predict_Robjs_rna/out_integration_ciliary_combined/obj_integ_ciliary_combined_rs_0.1_ds_30_ciliary_combined.rds")
prefix <- "ciliary_combined_degs"
dirOut <- "/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scRNA_HAS/scRNA_DEGs"

objA.markers <- runDEGspng(objectA, prefix,dirOut)


