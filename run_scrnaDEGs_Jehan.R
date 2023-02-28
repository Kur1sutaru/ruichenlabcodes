# title: draw vlnplot with p-values using ggpubr package


### Function: Vlnplot with p-values
drawVlnplotwPval<- function(gene_signature, file_name, test_sign1,test_sign2,test_sign3,test_sign4, plabel="p.signif"){
  plot_case1 <- function(signature, y_max = NULL)
  {
    VlnPlot(seurat, features = signature,
            pt.size = 0.1,group.by = "sample", y.max = y_max) + 
      stat_compare_means(comparisons = test_sign1, label = plabel,label.y = (y_max_list[[gene]] ))+ labs(subtitle=paste0("cluster=",setClust))+
      stat_compare_means(comparisons = test_sign2, label = plabel,label.y = (y_max_list[[gene]] + 0.5))+
      stat_compare_means(comparisons = test_sign3, label = plabel,label.y = (y_max_list[[gene]] + 1))+
      stat_compare_means(comparisons = test_sign4, label = plabel,label.y = (y_max_list[[gene]] + 1.5))
  }
  
  plot_list  <- list()
  y_max_list <- list()
  for (gene in gene_signature) 
  { 
    yvalue<- max(FetchData(seurat, vars = gene))
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]])
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 4) )
  }
  print(y_max_list)
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, ".png")
  ggsave(file_name, width = 14, height = 8)
}
###================================================


#library
library(Seurat)
library(ggplot2)
library(ggpubr)


# user's parameters
dirIn  <- "/path/to/input_rds/"
finrds <- "obj_DS_merged_1.rds"
gene_sig <- c("Il17a","Dcn") # get target genes
fileName <- "Vlnplot_comparsion_DS_genes" # outfilename
setClust <- "1" # set cluster number to check: e.g between 1-19


# reorder samples in Vlnplot
setwd(dirIn)
seurat<- readRDS(finrds)
seurat$sample <- factor(x = seurat$sample, levels = c('NS', 'DS1','DS3','DS5','DS10'))
print(table(seurat$seurat_clusters))

# subset
if(exists("setClust")){
  Idents(seurat)<-seurat$seurat_clusters
  seurat<- subset(seurat, subset=seurat_clusters==setClust)
}

# set comparison
comparisons1 <- list(c("NS", "DS1")) # list(names(table(objA$sample))) 
comparisons2 <- list(c("NS", "DS3")) 
comparisons3 <- list(c("NS", "DS5")) 
comparisons4 <- list(c("NS", "DS10")) 

setwd(dirOut)

# draw vlnplot with p-value from wilcox.test
# options: 
# plabel=("p.signif","p.format","p.adj",)
if(exists("setClust")){fileName=paste0(fileName,"_cl_",setClust)}
drawVlnplotwPval(gene_signature = gene_sig, file_name = fileName, test_sign1 = comparisons1, test_sign2 =comparisons2, test_sign3 =comparisons3, test_sign4 =comparisons4,plabel="p.format") 







### DEG analysis by subset of celltype
# test DEGs

dirIn  <- "/input/rds/path/" # input data(.rds) path
dirOut <- "out_DEGs" #output dir name
fin    <- "GDT.rds" # R obj file name 
prefix <- "gdT" # data set name


# library
library(Seurat)
library(ggrepel)

# function
setDir <- function(dirIn, dirOut){
  setwd(dirIn)
  dirIn <- getwd()
  
  if(!file.exists(dirOut))
    dir.create(dirOut)
  setwd(dirOut) 
  dirOut <- getwd()
  cat("dirOut: ",dirOut,"\n")
  return(dirOut)
}

# dirOut
dirOut<-setDir(dirIn,dirOut)

# load data
setwd(dirIn)
objA <- readRDS(fin)
objA$celltype.subset <- paste(objA$seurat_clusters, objA$orig.ident, sep = "_")
table(objA$celltype.subset)

# DEGs
DefaultAssay(object = objA) <- "RNA"
Idents(objA) <- objA$celltype.subset
objA.marker  <- FindMarkers(objA, ident.1 = "2_Pinkie",  ident.2 = "2_B6", 
                            logfc.threshold = 0.25,min.pct = 0.1, test.use = "wilcox") 
# option: wilcox and others https://satijalab.org/seurat/articles/de_vignette.html

head(sub1.marker)
objA.marker$Gene <- rownames(objA.marker)
dim(objA.marker)

# save DEG output
setwd(dirOut)
write.csv(objA.marker,paste0("DEG_",prefix,"_bysubgroup.vs.cluster.csv"), row.names = T)


# draw volc plot for significant DEGs
# cut-off 
adjpValue <- 0.001
fcvalue   <- 0.5

# vol plot
DEGout    <- objA.marker
sTitle    <- prefix
DEGout$Significant <- ifelse(DEGout$p_val_adj < adjpValue &  abs(DEGout$avg_log2FC) > fcvalue, "FDR < 0.001", "Not Sig")
head(DEGout)

p1 <- ggplot(DEGout, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant)) + ggtitle(paste0("Volcano plot of DEGs:",sTitle)) +
  scale_color_manual(values = c("red", "black")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +geom_text_repel(
    data = DEGout %>% filter (abs(DEGout$avg_log2FC) > fcvalue & DEGout$p_val_adj < adjpValue ),
    aes(label = Gene),size = 5, box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))+
  geom_vline(xintercept=c(-1*fcvalue, fcvalue), col="red") 

# save plot
pdf(paste0("plot_volcano_DEG_",prefix,".pdf"))
  print(p1)
dev.off()











# run DEGs
objA <- runDEGsubsetfrom2CombinedData(objA,prefix,dirOut = dirOut,optAssay = "RNA")


runDEGsubsetfrom2CombinedData <- function(objectA,prefix,dirOut,optAssay){
  # set dir
  subDir <- paste0("out_DEGs_bySubset_",prefix)
  setwd(dirOut)
  dir.create(subDir)
  setwd(subDir)

  
  Idents(objectA) <- objectA$celltype.subset
  
  # get celltype
  print(table(objectA$celltype.subset))
  ctypelist <- names(table(objectA$seurat_clusters))
  ctypelist
  clist <- names(table(objectA$celltype.subset))
  outDegs <- data.frame("Celltype","set1","set2","#DEGs")
  for (i in 1:length(ctypelist)) {
    print(ctypelist[i])
    cid   <- ctypelist[i]
    cpair <- grep(cid,clist,value = T)
    print(cpair)
    cpair.marker <- FindMarkers(objectA, ident.1 = cpair[1], ident.2 = cpair[2], verbose = FALSE, assay=optAssay)
    cat("cpair: ",cpair,": DEG#: ",dim(cpair.marker)[1],"\n")  
    write.csv(cpair.marker,paste0("DEGs_",cpair[1],".vs.",cpair[2],".csv"))
    eachdeg <- c(ctypelist[i],cpair[1],cpair[2],dim(cpair.marker)[1])
    outDegs <- rbind(outDegs,eachdeg)
    rm(eachdeg)
    saveRDS(cpair.marker,paste0("DEGs.marker_",prefix,"_",cid,".rds"))
  }
  colnames(outDegs) <- outDegs[1,]
  outDegs <- outDegs[-1,]
  write.csv(outDegs,paste0("summary_DEGs_",prefix,"_bySubset.celltype.csv"),row.names = F)
  
  return(objectA)
}

# test DEGsubset
# integ.p90 <- runDEGsubsetfrom2CombinedData(integ.p90,prefix = paste0(prefix,"_test"),dirOut = dirOut)















