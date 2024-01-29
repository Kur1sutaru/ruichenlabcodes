#### test scATAC data QC ====
# library
library(Signac)
library(EnsDb.Mmusculus.v79)
library(Seurat)
library(ggplot2)
library(GenomeInfoDb)
library(cowplot)
library(rmarkdown)
library(knitr)

sDir="/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/"
source(paste0(sDir,"R_functions_scATAC.R"))
set.seed(1234)

# user's parameters
#dirIn <- "/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scATAC_HAS/22_0688_Cornea_ATAC/outs/"
#dirOut <- "out_QC_signac"
#prefix <- "22_0688_Cornea"

# get parameters
args = commandArgs(trailingOnly=TRUE)
print(args)

if (length(args)<=1) {
#  stop("please type 3 parameters:dirIn, prefix(sample name), and genomeVersion('hg38', or 'mm10')", call.=FALSE)
  stop("please type 2 parameters:dirIn, prefix(sample name)", call.=FALSE)
} else {
  # default output file
  dirIn   <- args[1] # input folder names
  prefix  <- args[2] # sample names
  #genomeV <- args[3] # genomeversion- "hg38", or "mm10"
  #dirOut <- args[3] # path to create output folder names
  cat("dirIn: ",dirIn,"\n")
  cat("prefix: ",prefix,"\n")
  #cat("genomeV: ",genomeV,"\n")
}


# create output dir
dirOut <-paste0("out_QC_signac_",prefix)
#dirOut <- setDir(dirIn,dirOut)
dirOut <- setDir(dirIn=getwd(),dirOut)
genomeV <-"mm10"


#if(F){
# 1) generate signac object
objAt <- runcreateSignacObj(dirIn,dirOut,prefix,genomeV)

# 2) run QC
objAt<-runQCSignacObj(objAt,dirOut,prefix)

# 3) cell proportion
runCellProportion(objectA=objAt,prefix,dirOut,cType="ID")

# 3.2) feature plot
#objAt<-readRDS("/storage/chenlab/Data/HCA/HAS2/src/atacQC/out_QC_signac_ATAC-22-0183-RPE-WOLF/obj_HAS_ATAC-22-0183-RPE-WOLF_geneActivity_pf1000_11977.27_pct15_br0.05_ns4_te2_n6875.rds")
#dirOut <-"/storage/chenlab/Data/HCA/HAS2/src/atacQC/out_QC_signac_ATAC-22-0183-RPE-WOLF/"
DefaultAssay(objAt) <- "RNA"
checkTSSfeatureHS(objA=objAt,prefix,dirOut)

#}
# 4) generate report
# parameters for rendering=report_form.Rmd
callRmd <- paste0(sDir,"/report_form_QC.Rmd")
render_report(dirIn,prefix,dirOut,callRmd)
#}

cat("The SC process was done! : ", prefix,"\n")

