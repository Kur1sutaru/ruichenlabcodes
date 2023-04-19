library(Seurat)
library(dplyr)
library(ggplot2)
source("/storage/chen/home/sangbaek/pipeline/pipeline_taco/R_script/Retina_function_cellmarker.R",encoding="UTF-8")

# user's parameter
dirIn  <-"/path/to/rds_file"
finrds <-"obj_name.rds"
dirOut <-"out_assignCelltype"
prefix <-"sampleName"
newCellTypes <- c("Rod",	"Rod",	"Cone",	"BC.on",	"MG",	"RBC",	"AC.HC.RGC",	"BC.off",	"Unknown",	"AC.HC.RGC")

# create out dir
dirOut <-setDir(dirIn,dirOut)

# Assign cell types
obj_data <- runAssignCelltype(objectA = obj_data,prefix = prefix,newCellTypes = newCellTypes)


