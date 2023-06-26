# title: run_DARscATAC_v2.R
# parameters: finrds, prefix, setGrp
# usage: Rscript --vanilla run_DARscATAC_v2.R /path/atac.rds prefix setGrp(e    .g.predicted.id)
 if(F){
   finrds <- "/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/data/test_merge_atac_objs/out_merged_atac_objs_atac_TM2/obj_annotAtac_atac_TM2_wcompeaks_wharmony.rds"
   prefix <- "TM_combined"
   setGrp <- "predicted.id" # choose a group name(category col name) in the meta data to perform DAR
   mpct   <- 0.2 # default: 0.1
 }

# check parameters
 args = commandArgs(trailingOnly=TRUE)
 print(args)
 if (length(args)!=3) {
   stop("please type 3 parameters:input rds file(w/path), prefix, category name to calculate DAR
        (e.g. seurat_clusters, predicted.id, or celltype) ", call.=FALSE)
 } else {

# default output file
   finrds  <- args[1]
   prefix  <- args[2]
   setGrp  <- args[3]
 }

# mcpt: min of pct per feature
 mpct   <- 0.2 # default: 0.1


# library
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
set.seed(1234)

# get sys.time: format(Sys.time(), "%d_%b_%Y_t%H.%M")

# set src dir
 srcDir <-"/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/"
 source(paste0(srcDir,'R_functions_scATAC.R'))

# load data
 atac <- readRDS(finrds)
 print(atac)

# create out dir
 dirIn  <- getwd()
 dirOut <- paste0("out_atac_DAR_",prefix)
 dirOut <- setDir(dirIn, dirOut)

# check data
 # plot
 #p1<-DimPlot(atac,group.by="seurat_clusters",label = TRUE, pt.size = 0.1) + NoLegend()
 #p2<-DimPlot(pbmc,group.by="celltype",label = TRUE, pt.size = 0.1) #scpred_prediction
 p2<-DimPlot(atac,group.by="predicted.id",label = TRUE, pt.size = 0.1) #scpred_prediction
 #p3<-DimPlot(pbmc,group.by="seqType",label = F, pt.size = 0.1) + labs(caption = paste0("- Data: ",prefix))
 #p3<-DimPlot(atac,group.by=Idents(atac),label = F, pt.size = 0.1) + labs(caption = paste0("- Data: ",prefix))

setwd(dirOut)
 #pdf(paste0(prefix, "_2.2_featureplot_celltype.markers.pdf"),width = 9, height =15 )
 png(paste0("01_umap_",prefix,"_2.png"), width=900, height = 800, units="px")
 print(p2) #/p3)
 dev.off()

# call peaks
# coverage plot
 setwd(dirOut)
 png(paste0("02_1_plot_coverage_",prefix,".png"), width=500, height = 400, units="px")
 CoveragePlot(object = atac, group.by = "predicted.id", region = "Ccl2")
 dev.off()



