# check cell type  based on known markers for HAS data
# usage: Rscript --vanilla /storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/check_marksersHAS6tissue_v2.R /storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scRNA_HAS/integration/predict_Robjs_rna/out_integration_ciliary_combined/obj_integ_ciliary_combined_rs_0.1_ds_30_ciliary_combined.rds ciliary_combined


# library
source("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/R_functions_scATAC.R")
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(reticulate)

#user parameters
#objectA   <- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scRNA_HAS/integration/predict_Robjs_rna/out_integration_ciliary_combined/obj_integ_ciliary_combined_rs_0.1_ds_30_ciliary_combined.rds")
#prefix    <- "ciliary_combined"

# check parameters
args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)!=2) {
  stop("please type 2 parameters:input file(w/path)", call.=FALSE)
} else {
  # default output file
  finrds  <- args[1]
  prefix  <- args[2]
}

objectA <- readRDS(finrds)
# load marker Sets
markerSet <- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/marker_genes_HAS/obj_marker_HAS_list_6T.rds")
dirOut    <- getwd()

# check markers for 6 tissue types  
setwd(dirOut)
checkMarkersHsHas(objectA,markerSet,prefix,dirOut)

run_dotplot4HAScelltypeMkers(objectA,prefix)

print("Process was done!")

