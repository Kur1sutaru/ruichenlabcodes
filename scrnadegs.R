### set dir -----
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

### DEGs ----
runDEGspng <- function(objectA, prefix, dirOut){
# set dir
subDir <- paste0("out_DEGs_",prefix)
setwd(dirOut)
dir.create(subDir)
setwd(subDir)

# DEGs
#Identify conserved cell type markers: for two groups
DefaultAssay(object = objectA) <- "RNA"
# scale data
objectA <- ScaleData(object = objectA)

#find differentially expressed gene markers for every cluster compared to all remaining cells
objA.markers <- FindAllMarkers(object = objectA, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
head(objA.markers)

# select top 20, 50 genes of each cluster
aTop2<- objA.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
aTop5<- objA.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
aTop10 <- objA.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
aTop50 <- objA.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
aTop100 <- objA.markers %>% group_by(cluster) %>% top_n(100, avg_log2FC)

# write table - top50 markers
write.csv(aTop50, paste0("out_top50markers_eachCluster_",prefix,".csv"),row.names= F)
write.csv(aTop100, paste0("out_top100markers_eachCluster_",prefix,".csv"), row.names= F)

# save - DEG.obj
saveRDS(objA.markers,paste0("objA.marker_",prefix,".rds"))


# Heatmap
p1<-DoHeatmap(object = objectA, features = aTop10$gene) + NoLegend()
#pdf(paste0("1_1_plot_heatmap_top10_",prefix,".pdf"), height = 20)
png(paste0("1_1_plot_heatmap_top10_",prefix,".png"))
print(p1)
dev.off()

p2<-DoHeatmap(object = subset(objectA,downsample=500), features = aTop10$gene) + NoLegend()
#pdf(paste0("1.2_plot_heatmap_top10_by100cells_",prefix,".pdf"), height = 20)
png(paste0("1_2_plot_heatmap_top10_by100cells_",prefix,".png"))
print(p2)
dev.off()


 # plot_dot_top5
p3<- DotPlot(object = objectA, features = rev(x = unique(aTop5$gene)), dot.scale = 8, assay = "RNA",
 col.min = 0) + RotatedAxis()
#pdf(paste0("1.3.plot_dot_markers_top5_",prefix,"_1.pdf"), height=4,width = 15)
png(paste0("1_3_plot_dot_markers_top5_",prefix,".png"), height=400,width = 1200)
print(p3)
dev.off()

p4<- DotPlot(object = objectA, features = rev(x = unique(aTop10$gene)), dot.scale = 8, assay = "RNA",
 col.min = 0) + RotatedAxis()
#pdf(paste0("1.4.plot_dot_markers_top10_",prefix,".pdf"), height=4,width = 25)
png(paste0("1_4_plot_dot_markers_top10_",prefix,".png"), height=400,width = 2000)
print(p4)
dev.off()

return(objA.markers)
}

