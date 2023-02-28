
# input data-DEGs per celltype
# done
if(F){
  dirIn <- "/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/"
  prefix<- "AC"
  dirIn <- "/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/"
  prefix<- "BC"
  dirIn <- "/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/"
  dirIn <- "/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/"
  dirIn <- "/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/"
  dirIn <- "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/"
}
dirIn <- "/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/"


dirOut <- "/storage/singlecell/sangbaek/307_10xscRNA_RPE/09_compare_RetinaData/"
fin1 <- "Early_compare_mod_region.csv"
fin2 <- "Early_fit_coefs_region_models.csv"

# library
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(reticulate)

# set src dir
srcDir <-"/storage/chen/data_sharefolder/22_10x_hs_AnteriorSegment_data/scAtacQC/"
source(paste0(srcDir,'R_functions_scATAC.R'))


# get common gList
setwd(dirIn)
commonGListFovPer <- doIntersectDEGsRetinaCell(fin1, fin2, prefix, dirOut)


# library
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(reticulate)

# set src dir
srcDir <-"/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/"
source(paste0(srcDir,'R_functions_scATAC.R'))

# load data-DEGs
# Function-celltype DEGs in a celltype
# usage: commonGListFovPer <- doIntersectDEGsRetinaCell(fin1, fin2, prefix, dirOut)
setwd(dirIn)
doIntersectDEGsRetinaCell<-function(fin1, fin2, prefix,dirOut){
  df <- read.csv(fin1)
  df <- df[complete.cases(df),]
  df <- df[df$q_value < 0.01,]
  
  df2 <- read.csv(fin2)
  df2 <- df2[complete.cases(df2),]
  df2 <- df2[(df2$q_value < 0.01)&(df2$term=="RegionPeripheral")&(df2$status=='OK')&(abs(df2$normalized_effect)>1),]
  print(dim(df2))

  # common DEGs-w/o region vs w/ region
  gsdf<- intersect(df$gene_short_name,df2$gene_short_name)
  gsdf2<-merge(df,df2,by.x = "gene_short_name", by.y = "gene_short_name")
  ncg <-dim(gsdf2)[1]
  setwd(dirOut)
  write.csv(gsdf2, paste0("DEGs_common_foveal_vs_peri_retina_",prefix,"_n_",ncg,".csv"))
  return(gsdf2)
}



setwd(dirIn)
df <- read.csv(fin1)
df <- df[complete.cases(df),]
df <- df[df$q_value < 0.01,]

df2 <- read.csv(fin2)
df2 <- df2[complete.cases(df2),]
df2 <- df2[(df2$q_value < 0.01)&(df2$term=="RegionPeripheral")&(df2$status=='OK')&(df2$normalized_effect>1),]


# common DEGs-w/o region vs w/ region
gsdf<- intersect(df$gene_short_name,df2$gene_short_name)
gsdf2<-merge(df,df2,by.x = "gene_short_name", by.y = "gene_short_name")
write.csv(gsdf2, paste0("DEGs_common_foveal_vs_peri_retina_",prefix,".csv"))


# print gene list
# for (x in (intersect(df$gene_short_name,df2$gene_short_name))){cat(x,"\n")}