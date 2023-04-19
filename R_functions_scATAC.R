
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


## function - runQCscATAC.R
# parameter: dirIn-dir for .bam, dirOut-output folder, prefix
# QC parameters
# minPRK=3000  # minimum of fragments in peaks
# maxPRK=20000 # max of fragments in peaks
# PctRiP=15  # Fraction of fragments in peaks    
# blr=0.05 # Ratio reads in genomic blacklist regions
# nsr=4  # 
# tssE=2 # Transcriptional start site (TSS) enrichment score

# usage: objAt <- runQCscATAC(dirIn, dirOut,prefix,genomeV)
runcreateSignacObj <- function(dirIn,dirOut,prefix,genomeV){
  #library(Signac);library(EnsDb.Hsapiens.v86);library(EnsDb.Mmusculus.v79)
  # QC parameters
  #minPRK=3000  # minimum of fragments in peaks
  #maxPRK=20000 # max of fragments in peaks
  #PctRiP=15  # Fraction of fragments in peaks    
  #blr=0.05 # Ratio reads in genomic blacklist regions
  nsr=4  # 
  tssE=2 # Transcriptional start site (TSS) enrichment score
  #genomeV="hg38" #options: "hg38","mm10","rn6" and others
  
  # set dir
  if(is.null(dirOut)){dirOut <- "out_scATAC_signac"}
  dirOut <- setDir(dirIn=dirIn, dirOut = dirOut)
  
  # input data
  setwd(dirIn)
  counts   <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")
  metadata <- read.csv(file = "singlecell.csv", header = TRUE,  row.names = 1)
  
  chrom_assay <- CreateChromatinAssay(counts = counts,  sep = c(":", "-"),genome = genomeV,  
                                      fragments = 'fragments.tsv.gz', min.cells = 1,  min.features = 10)
  
  objAt <- CreateSeuratObject(counts = chrom_assay,  assay = "peaks",  meta.data = metadata)
  objAt
  objAt$orig.ident <- prefix
  
  # annotation
  # extract gene annotations from EnsDb
  if(genomeV=="hg38"){
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)} else if(genomeV=="mm10"){
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)} else {
        stop("check genome!")}
  
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  # add the gene information to the object
  Annotation(objAt) <- annotations
  
  
  # compute nucleosome signal score per cell
  objAt <- NucleosomeSignal(object = objAt)
  
  # compute TSS enrichment score per cell
  objAt <- TSSEnrichment(object = objAt, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  objAt$pct_reads_in_peaks <- objAt$peak_region_fragments / objAt$passed_filters * 100
  objAt$blacklist_ratio <- objAt$blacklist_region_fragments / objAt$peak_region_fragments
  
  objAt$high.tss <- ifelse(objAt$TSS.enrichment > tssE, 'High', 'Low')
  p1<-TSSPlot(objAt, group.by = 'high.tss') + NoLegend() + ggtitle(paste0("TSS enrichment: ",prefix))+labs(subtitle=prefix)
  
  objAt$nucleosome_group <- ifelse(objAt$nucleosome_signal > nsr, paste0('NS > ',nsr), paste0('NS < ',nsr))
  p2<- FragmentHistogram(object = objAt, group.by = 'nucleosome_group') +labs(subtitle=prefix)
  p3<-VlnPlot(object = objAt,
              features = c('pct_reads_in_peaks', 'peak_region_fragments',
                           'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),pt.size = 0,  ncol = 5)  +labs(subtitle=prefix)
  
  p4 <- plot_grid(p1,p2,ncol=1)
  
  # save RDS
  setwd(dirOut)
  ncells <-dim(objAt)[2]
  saveRDS(objAt, paste0("obj_HAS_",prefix,"_unfilt_n",ncells,".rds"))
  
  # plot-tss
  setwd(dirOut)
  #pdf(paste0("01_1_plot_tss_",prefix,".pdf"))
  png(paste0("01_1_plot_tss_",prefix,".png"))
  print(p4)
  dev.off()
  
  setwd(dirOut)
  png(paste0("01_2_plot_QC_",prefix,"_unfilted.png"),height = 300, width = 1000)
  print(p3)
  dev.off()
  
  return(objAt)
}  
  



# function-qc summary table from scATAC Object
# usage: qcsumTable <- runACsumT(objAt,prefix)

runQCsumT<-function(objAt,prefix){
  prp<- objAt$pct_reads_in_peaks
  sd.prp <- sd(prp) *3 + mean(prp)
  prps <- summary(prp)
  
  pf <- objAt$peak_region_fragments
  sd.pf <- sd(pf) *3 + mean(pf)
  #sd.pf2 <- sd(pf) * -3 + mean(pf)
  pfout <- sum(pf>sd.pf)
  pfs <- summary(pf)
  
  te <-objAt$TSS.enrichment
  sd.te <- sd(te) *3 + mean(te)
  tes <-summary(te)
  
  nss <-objAt$nucleosome_signal
  sd.nss <- sd(nss)*3 + mean(nss)
  nsss<- summary(nss)
  
  sumAt <-round(rbind(prps,pfs,tes,nsss),2)
  rName <-c("pct_reads_in_peaks","peak_region_fragments","TSS.enrichment","nucleosome_signal")
  rownames(sumAt)<-rName
  sd3cutoff=round(c(sd.prp, sd.pf, sd.te, sd.nss),2)
  sumAt<-cbind(sumAt,sd3cutoff)
  write.csv(sumAt,paste0("summary_QC_statistics_unFilt_",prefix,".csv"))
  print(sumAt)
  return(sumAt)
}

## function - run QC for scATAc data
# minPRK=3000  # minimum of fragments in peaks
# maxPRK=20000 # max of fragments in peaks
# PctRiP=15  # Fraction of fragments in peaks    
# blr=0.05 # Ratio reads in genomic blacklist regions
# nsr=4  # 
# tssE=2 # Transcriptional start site (TSS) enrichment score


# cut-off parameter
runQCSignacObj <- function(objAt,dirOut,prefix){
  # default parameters
  minPRK =1000
  PctRiP=15 
  blr=0.05 
  nsr=4  
  tssE=2 
  
  # QC summary table
  
  qcsumTable <- data.frame(runQCsumT(objAt,prefix))
  maxPRK <-qcsumTable$sd3cutoff[2]
  if(maxPRK <10000)(maxPRK=10000)
  cat("maxPRK= ",maxPRK,"\n")
  
  # filtering
  objAt <- subset(
    x = objAt,
    subset = peak_region_fragments > minPRK &
      peak_region_fragments < maxPRK &
      pct_reads_in_peaks > PctRiP &
      blacklist_ratio < blr &
      nucleosome_signal < nsr &
      TSS.enrichment > tssE
  )
  qcCells<-dim(objAt)[2]
  p4<-VlnPlot(object = objAt,
              features = c('pct_reads_in_peaks', 'peak_region_fragments',
                           'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
              pt.size = 0,  ncol = 5) + labs(subtitle=prefix)
  setwd(dirOut)
  png(paste0("01_3_plot_QC_",prefix,"_filted_n",qcCells,".png"),height = 300, width = 1000)
  print(p4)
  dev.off()
  
  print(objAt)
  # save RDS
  setwd(dirOut)
  ncells <-dim(objAt)[2]
  saveRDS(objAt, paste0("obj_HAS_",prefix,"_qc_pf",minPRK,"_",maxPRK,"_pct",PctRiP,"_br",blr,"_ns",nsr,"_te",tssE,"_n",ncells,".rds"))
  
  # normalization, linear dimensional reduction
  objAt <- RunTFIDF(objAt)
  objAt <- FindTopFeatures(objAt, min.cutoff = 'q0')
  objAt <- RunSVD(objAt)
  
  p5<-DepthCor(objAt)+ labs(subtitle=prefix)
  setwd(dirOut)
  #pdf(paste0("02.1_plot_depth_",prefix,".pdf"),width = 6, height = 4)
  png(paste0("02_1_plot_depth_",prefix,".png"),width = 600, height = 500)
  print(p5)
  dev.off()
  
  # Non-linear dimension reduction and clustering => removed lsi 1 dim
  objAt <- RunUMAP(object = objAt, reduction = 'lsi', dims = 2:30)
  objAt <- FindNeighbors(object = objAt, reduction = 'lsi', dims = 2:30)
  objAt <- FindClusters(object = objAt, verbose = FALSE, algorithm = 3)
  p6 <- DimPlot(object = objAt, label = TRUE) + NoLegend()+ labs(subtitle=prefix)
  #pdf(paste0("02.2_plot_umap_",prefix,".pdf"),width = 6, height = 4)
  png(paste0("02_2_plot_umap_",prefix,".png"),width = 600, height = 500)
  print(p6)
  dev.off()
  
  # Create a gene activity matrix
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  gene.activities <- GeneActivity(objAt)
  
  objAt[['RNA']] <- CreateAssayObject(counts = gene.activities)
  objAt <- NormalizeData(
    object = objAt,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(objAt$nCount_RNA)
  )
  
  objAt <- RunUMAP(object = objAt, reduction = 'lsi', dims = 2:30)
  objAt <- FindNeighbors(object = objAt, reduction = 'lsi', dims = 2:30)
  objAt <- FindClusters(object = objAt, verbose = FALSE, algorithm = 3)
  p7 <- DimPlot(object = objAt, label = TRUE) + NoLegend()+ labs(subtitle=prefix)
  #pdf(paste0("02.3_plot_umap_afterNorm_geneActivity_",prefix,".pdf"),width = 5, height = 4)
  png(paste0("02_3_plot_umap_afterNorm_geneActivity_",prefix,".png"),width = 600, height = 500)
  print(p7)
  dev.off()
  
  # save RDS
  setwd(dirOut)
  ncells <-dim(objAt)[2]
  saveRDS(objAt, paste0("obj_HAS_",prefix,"_geneActivity_pf",minPRK,"_",maxPRK,
                        "_pct",PctRiP,"_br",blr,"_ns",nsr,"_te",tssE,"_n",ncells,".rds"))
  
  
  return(objAt)
  cat("Process was done!- runQCSignacObj","\n")
  
} # end; runQCscATAC

### cell proportion ----
cType <- "ID" #option: ID,Clust,subset 
runCellProportion <- function(objectA,prefix,dirOut,cType){
  # cell proportion
  if(cType=="ID"){
    metaTbl  <- data.frame(table(Idents(objectA)))
  }else if(cType=="Clust"){
    metaTbl  <- data.frame(table(objectA$seurat_clusters))
  }else if(cType=="subset"){
    metaTbl  <- data.frame(table(objectA$subset))
  }else {print("choose cluster type:(ID or Clust")}
  prop2pr  <- data.frame(proportion=round(prop.table(metaTbl$Freq)*100,1))
  metaTblA <- cbind(metaTbl,prop2pr)
  colnames(metaTblA) <- c('cluster','cells','Proportion')
  
  setwd(dirOut)
  write.csv(metaTblA, paste0(prefix,"_proportion_clusters.csv"),row.names= F)
  
  # cell table
  if(!(sum(names(objectA@meta.data) == "Cluster"))){
    clustering.table <- table(Idents(objectA), objectA$seurat_clusters)
    head(clustering.table)
    setwd(dirOut)
    write.csv(clustering.table, paste0(prefix,"_summary_table_clusters_celltype_",prefix,".csv"),row.names= T)
  }
   
  # plot-bar
  #pdf(paste0(prefix,"_plot_Barplot_cell.proportion.pdf"),width = 20, height = 10)
  png(paste0(prefix,"_plot_Barplot_cell.proportion.png"),width = 1200, height = 800)
  print(ggplot(metaTblA, aes(x=cluster, y=Proportion, fill=cluster)) + ylab("Proportion (%)") +xlab("Cluster")+
          geom_bar(stat="identity") + geom_text(aes(label=cells), vjust=-0.3, size=15) +
          scale_x_discrete(limits = rev(levels(metaTblA$cluster)))+ coord_flip()+
          theme(text = element_text(size=20)))
  dev.off()
  cat("Process was done! \n")
}
# test run: cell proportion
# runCellProportion(obj.WTp14,"test",dirOut = dirOut,cType = "ID")

### function : generare report
#render_report = function(dirIn, prefix, callRmd,dirOut) {
render_report = function(dirIn,prefix,dirOut,callRmd) {
  rmarkdown::render(
    callRmd, params = list(
#      dirIn = dirIn,
      prefix=prefix,
#      finsum = finsum,
      dirOut=dirOut),
    output_file = paste0(dirOut,"/Report_QC_", prefix, ".html"),"html_document")
}


# Function: runCoembedscRnaAtac - do co-embedding for celltype prediction
# require: library(EnsDb.Hsapiens.v86)
# rnaCellType="celltype" # celltype info column name in scRNA
# parameters: fragPath-frgment file path
runCoembedscRnaAtacHg38 <- function(scRNA,scAtac,dirOut,prefix,rnaCellType="celltype",fragPath){
  print(scRNA); print(scAtac)
  scRNA$seqType  <- "scRNA"
  scAtac$seqType <- "scATAC"
  
  # update fragment file path
  frags <- Fragments(scAtac)  # get list of fragment objects
  Fragments(scAtac) <- NULL  # remove fragment information from assay
  newpath <- paste0(fragPath,"/fragments.tsv.gz")
  frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
  Fragments(scAtac) <- frags

  # Perform standard analysis of each modality independently RNA analysis
  DefaultAssay(object = scRNA) <- "RNA"
  scRNA <- NormalizeData(scRNA)
  scRNA <- FindVariableFeatures(scRNA)
  scRNA <- ScaleData(scRNA)
  scRNA <- RunPCA(scRNA)
  scRNA <- RunUMAP(scRNA, dims = 1:30)
  
  # ATAC analysis add gene annotation information
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  Annotation(scAtac)  <- annotations
  
  
  # We exclude the first dimension as this is typically correlated with sequencing depth
  scAtac <- RunTFIDF(scAtac)
  scAtac <- FindTopFeatures(scAtac, min.cutoff = "q0")
  scAtac <- RunSVD(scAtac)
  scAtac <- RunUMAP(scAtac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  # plot
  p1 <- DimPlot(scRNA, group.by = rnaCellType, label = TRUE) + NoLegend() + ggtitle("RNA")
  p2 <- DimPlot(scAtac, group.by = "seqType", label = FALSE) + NoLegend() + ggtitle("ATAC")
  setwd(dirOut)
  png(paste0("1_plot_umap_RNA_ATAC_",prefix,".png"),width=800, height=300)
    print(p1 + p2)
  dev.off()
  
  # Identifying anchors between scRNA-seq and scATAC-seq datasets
  # quantify gene activity
 # setwd(fragPath)
  gene.activities <- GeneActivity(scAtac, features = VariableFeatures(scRNA))
  
  # add gene activities as a new assay
  scAtac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
  
  # normalize gene activities
  DefaultAssay(scAtac) <- "ACTIVITY"
  scAtac <- NormalizeData(scAtac)
  scAtac <- ScaleData(scAtac, features = rownames(scAtac))
  
  # Identify anchors
  transfer.anchors <- FindTransferAnchors(reference = scRNA, query = scAtac, features = VariableFeatures(object = scRNA),
      reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
  
  # Annotate scATAC-seq cells via label transfer
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = scRNA$celltype,
      weight.reduction = scAtac[["lsi"]], dims = 2:30)
  
  scAtac <- AddMetaData(scAtac, metadata = celltype.predictions)
  
  # save RDS
  setwd(dirOut)
  saveRDS(scAtac, paste0("obj_scATAC_wcelltypeTransfer_",prefix,".rds"))
  
  
  # plot
  #scAtac$annotation_correct <- scAtac$predicted.id == scAtac$seurat_annotations
  p3 <- DimPlot(scAtac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
  #p3.2 <- DimPlot(scAtac, group.by = rnaCellType, label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
  
  setwd(dirOut)
  png(paste0("2_plot_umap_RNA_ATAC_labeTransfer_",prefix,".png"))
    print(p3)
  dev.off()
  
  # Co-embedding scRNA-seq and scATAC-seq datasets
  genes.use <- VariableFeatures(scRNA)
  refdata   <- GetAssayData(scRNA, assay = "RNA", slot = "data")[genes.use, ]
  
  # imputation
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = scAtac[["lsi"]],dims = 2:30)
  scAtac[["RNA"]] <- imputation
  
  coembed <- merge(x = scRNA, y = scAtac)
  # add celltype
  #coembed$celltype2 <-c(as.character(scAtac$predicted.id),as.character(scRNA$celltype))
  
  # PCA and UMAP on this combined object to visualize the co-embedding 
  coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30)
  
  p4a <-DimPlot(coembed, group.by = "seqType")+ ggtitle("Co-embeded data")
  p4b <-DimPlot(coembed, group.by = "celltype")+ ggtitle("scRNA:celltype")
  p4c <-DimPlot(coembed, group.by = "predicted.id")+ ggtitle("scATAC:predicted celltype")
  
  
  #plot_annotation(title ='Co-embeded data, scRNA:celltype, scATAC:predicted celltype')
  png(paste0("3_plot_umap_RNA_ATAC_coembed_",prefix,".png"),width=1200,height=300)
    print(p4a+p4b+p4c)
  dev.off()
  
  # save RDS
  setwd(dirOut)
  saveRDS(coembed, paste0("obj_coembed_",prefix,".rds"))

  # cellproportion
  Idents(scAtac)<-scAtac$predicted.id
  figTitle=paste0("scATAC_",prefix)
  runCellProportion(objectA=scAtac,prefix,dirOut,cType = "ID")
  return(coembed)
}



## scRNA analysis-scPred

### Function : createscPredReference
# parameters: objectA(require 'celltype' info), dirOut,clssID(default:celltype) 
# usage: reference <- createscPredReference(objectA,dirOut,clssID, prefix)

createscPredReference <- function(objectA,dirOut,clssID,prefix){
  print(head(objectA@meta.data))
  #if(is.null(clssID)) {clssID <- "celltype"; print(clssID)} 
  
  # umap
  setwd(dirOut)
  p1<-DimPlot(object = objectA, reduction = "umap", group.by = clssID,label = T,label.size = 5)
  pdf(paste0("1_plot_umap_ref_",prefix,".pdf"),width = 8)
  print(p1)
  dev.off()
  
  # cell proportion
  cID <-which(names(objectA@meta.data)==clssID)
  head(objectA@meta.data[,cID])
  
  Idents(objectA) <- objectA@meta.data[,cID]
  cType <- "ID" #option: ID,Clust,subset,subset.all 
  runCellProportion(objectA=objectA,prefix=prefix,dirOut=dirOut,cType=cType)
  
  # train classifier
  objectA <- getFeatureSpace(objectA, pvar=clssID)
  objectA <- trainModel(objectA)
  get_probabilities(objectA) %>% head()
  probTable <- data.frame(get_probabilities(objectA))
  write.csv(probTable,paste0("summary_probTable_ref_",prefix,".csv"),row.names=T)
  
  # save RDS
  print(objectA)
  saveRDS(objectA, paste0("obj_reference_scPredTable_",prefix,".rds")) 
  
  sink(paste0("summary_scPred_Table_ref_",prefix,".out"))
  print(get_scpred(objectA)) #option: recompute_alignment = FALSE
  sink()
  
  p3<-plot_probabilities(objectA)
  png(paste0("2_plot_ref_probability_",prefix,".png"),width = 800,height = 800)
    print(p3)
  dev.off()

  
  return(objectA)
  
}

if(F){
  objectA<- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/data/test/out_scPred_ref_testSacilia/obj_reference_scPredTable_testSacilia.rds")
  prefix <- "test_cilia"
  dirOut<-"/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/data/test/out_scPred_ref_testSacilia/"
  list.file
  
}
### Function run scPred with ref and train data
# parameters: reference,query,prefix,dirOut, rsN=0.1,dn=20
# usage: query <- runscPred(objref=reference,objtrain=query,prefix,dirOut)

runscPred<-function(reference,query,prefix,dirOut){
  rsN=0.1
  dn=20
  # Training data- prediction =====
  setwd(dirOut)
  subDir <- paste0("out_cellType_prediction_",prefix)
  subDir <- setDir(dirIn=dirOut,dirOut=subDir)
  
  #lssvmRadial
  #query <- NormalizeData(query)
  query <- scPredict(query, reference)
  
  # save RDS
  setwd(subDir)
  saveRDS(query, paste0("obj_query_predicted_",prefix,".rds"))
  
  
  # umap- outputs
  query  <- RunUMAP(query, reduction = "scpred", dims = 1:dn)
  
  # cell clustering
  query <- FindNeighbors(object = query, reduction = "scpred", dims = 1:dn)
  query <- FindClusters(query, resolution = rsN)
  
  
  # cell proportion
  Idents(query) <- query$seurat_clusters
  cType <- "ID" #option: ID,Clust,subset,subset.all 
  runCellProportion(objectA=query,prefix=paste0(prefix,"_bycluster"),dirOut=subDir,cType=cType)
  
  p4 <-DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
  p5 <-DimPlot(query, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
  
  # plot-umap
  setwd(subDir)
  pdf( paste0("01_1_plot_umap_classf_scpred_",prefix,".pdf"),width = 10)
  print(p4) #| p2
  dev.off()
  
  pdf(paste0("01_2_plot_umap_classf2_compare_ref.vs.train_",prefix,".pdf"),width = 14)
  print(p5 |p4) #| p2
  dev.off()
  
  
  # save RDS
  setwd(subDir)
  saveRDS(query, paste0("obj_query_predicted_umap_",prefix,".rds"))
  
  #compare raw vs prediction
  outCross <- crossTab(query, "seurat_clusters", "scpred_prediction")
  write.csv(outCross,paste0("summary_coross_seuratcluster.vs.prediction_",prefix,".csv"))
  
  sink(file=paste0("summary_classifiers_",prefix,".out"))
  get_classifiers(reference)
  sink()
  
  return(query)
}



### Function: draw plots(dot,feature,vlnplot) fro marker set list
# parameters: objectA-Robj, markerSet-markers list format, prefix, dirOut
# usage:checkMarkersHsHas(objectA,markerSet,prefix,dirOut)

checkMarkersHsHas <- function(objectA,markerSet,prefix,dirOut){
  #markerSet #Mlist1
  markerSetNames <- names(markerSet)
  print(names(markerSet))
  # set dir
  subDir <- paste0("out_check_markerGeneExpression_",prefix)
  setwd(dirOut)
  dir.create(subDir)
  setwd(subDir)
  
  allGenes <- rownames(objectA)
  for (i in 1:length(markerSetNames)) {
    cat ("process_marker_set: ",i,"\n")
    mgene1   <- markerSet[[i]]
    mgeneIS1 <- toupper(mgene1)
    # set plot size
    iwid <- 8+length(mgeneIS1)*1/3
    if(length(mgeneIS1)/4 <=1){iwid2<-300*length(mgeneIS1)} else{iwid2= 1200}
    ihig2=300*ceiling((length(mgeneIS1)/4))
    ihig3=3*ceiling((length(mgeneIS1)/4))  
    # dot plot
    pdf(paste0(i,"_1_plot_dot_markers_",markerSetNames[i],"_",prefix,".pdf"),width = iwid)
    print(DotPlot(object = objectA,features = unique(mgeneIS1), dot.scale = 6, assay = "RNA",col.min = 0, scale.max = 80,
                  cols = c("gray","magenta1")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
            theme(panel.grid.major = element_line(colour = "grey95", size = (.1)))+ggtitle(paste0("Marker gene expression:",markerSetNames[i],"- ",prefix))) #+ coord_flip
    dev.off()
    
    # feature plot
    #pdf(paste0(prefix, "_2.2_featureplot_celltype.markers.pdf"),width = 9, height =15 )
    png(paste0(i,"_2_featureplot_celltype_markers_",markerSetNames[i],"_",prefix,"2.png"),width = iwid2, height = ihig2)
    print(FeaturePlot(object = objectA, features = mgeneIS1, pt.size=.1, min.cutoff = 0, ncol = 3)+
            labs(caption = paste0("TissueMarker: ",markerSetNames[i],"- Sample: ",prefix)))
    dev.off()
    # vlnplot
    #png(paste0(prefix,"_2_4_vlnplot_celltype.markers_",markerSetNames[i],"_",prefix,".png"),width = iwid2*1.5, height = ihig2))
    pdf(paste0(i,"_3_vlnplot_celltype_markers_",markerSetNames[i],"_",prefix,".pdf"),width = 12, height =ihig3)
    print(VlnPlot(object = objectA,features = mgeneIS1,pt.size = 0, ncol = 3)+
            labs(caption = paste0("TissueMarker: ",markerSetNames[i],"- Sample: ",prefix)))
    dev.off()
  }
  print(markerSet)
  print("Process was done!")
}



# function - getSubTissueMarker.R
# parameters: mkers-Robj with HAS list, tType-tissue or celltype e.g. TM
# usage: chkMkList2<- getSubTissueMarker(mkers,tType)

#tType='TM'
getSubTissueMarker<- function(mkers,tType){
  submkers<- dplyr::filter(mkers,Tissue %in% tType)[,c(2,3)]
  
  tNames2 <- names(table(submkers$celltype))
  chkMkList2 <- list()
  
  for (i in 1:length(tNames2)) {
    print(tNames2[i])
    eachtissue2 <- tNames2[i]
    submker2   <- dplyr::filter(submkers,celltype %in% eachtissue2)$marker
    chkMkList2[[i]] <- submker2
  }
  names(chkMkList2)<- tNames2
  print(head(chkMkList2))
  return(chkMkList2)
}


### function - run_dotplot4HAScelltypeMkers: check HAS markers showing celltype info on dot plot
# parameters: objectA-Robj with fullpath, prefix- sample_ID
# requirement: getSubTissueMarker function
run_dotplot4HAScelltypeMkers<- function(objectA,prefix){
  # load HAS marker object
  mkers    <- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/data/obj_mkers_HAS_mkerSets.rds")
  chkMkList<- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/data/obj_marker_HAS_list_6T.rds")
  
  # check tissue types
  tNames2 <- names(chkMkList)
  
  #setwd(dirOut)
  for (i in 1:length(tNames2))
  {
    tType <- tNames2[i]
    chkMkList2 <- getSubTissueMarker(mkers,tType)
    mgeneIS1   <- unique(unlist(chkMkList2))
    
    # dot plot for each celltype
    iwid <- (8+length(mgeneIS1)*1/3)*2
    pdf(paste0(i,"_1_plot_dot_markers_",tType,"_",prefix,".pdf"),width = iwid)
    print(DotPlot(object = objectA,features = chkMkList2, dot.scale = 6, assay = "RNA",col.min = 0, scale.max = 80,
                  cols = c("gray","magenta1")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
            labs(caption=paste0("-Data: ",prefix,", Tissue: ",tType)))#+
    #theme(panel.gratrid.major = element_line(colour = "grey95", size = (.1)))+ggtitle(paste0("Marker gene expression:",markerSetNames[i],"- ",prefix))) #+ coord_flip
    dev.off()
    #return(chkMkList)
  }
}

### run data integration ---------
runIntegrateObjects <- function(objList,dirOut,prefix,rsN,ds){
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = objList)
  
  # Perform integration
  integ.anchors   <- FindIntegrationAnchors(object.list = objList, anchor.features=features, dims = 1:ds)
  integ.combined  <- IntegrateData(anchorset = integ.anchors, dims = 1:ds)
  integ.combined
  
  DefaultAssay(integ.combined) <- "integrated"
  # Run the standard workflow for visualization and clustering
  integ.combined <- ScaleData(integ.combined, verbose = FALSE)
  integ.combined <- RunPCA(integ.combined, npcs = 30, verbose = FALSE)
  integ.combined <- RunUMAP(integ.combined, reduction = "pca", dims = 1:ds)
  integ.combined <- FindNeighbors(integ.combined, reduction = "pca", dims = 1:ds)
  integ.combined <- FindClusters(integ.combined, resolution = rsN)
  
  #integ.combined$subset <- gsub("MK.org.D173","D173",integ.combined$subset)
  integ.combined <- FindClusters(integ.combined, resolution = rsN)
  print(head(integ.combined@meta.data))
  
  # save rds
  setwd(dirOut)
  saveRDS(integ.combined,paste0("obj_integ_",prefix,"_rs_",rsN,"_ds_",ds,"_",prefix,".rds"))
  
  pdf(paste0("1_plot_umap_rs_",rsN,"_ds_",ds,"_",prefix,".pdf"))
  print(DimPlot(object = integ.combined, reduction = "umap",label=T,label.size = 8)) 
  print(DimPlot(object = integ.combined, reduction = "umap",group.by ="subset" ,label=F,label.size = 8))
  #print(DimPlot(object = integ.combined, reduction = "umap",group.by ="celltype" ,label=T,label.size = 6)) 
  #print(DimPlot(object = integ.combined, reduction = "umap",group.by ="celltype" ,label=F,label.size = 8)) 
  dev.off()
  
  nd <-dim(summary(objList))[1]
  nwid <- 300*nd
  png(paste0("2_plot_umap_rs_",rsN,"_ds_",ds,"_",prefix,".png"),width = nwid, height = 300)
   print(DimPlot(object = integ.combined, reduction = "umap",split.by ="subset" ,label=F,label.size = 5)) 
  dev.off()
  
  cat("Process was done!\n")
  return(integ.combined)
}


### run data integration ---------
runIntegrateObjectsATAC <- function(objList,dirOut,prefix,rsN,ds){
  
}




# Function - draw venndiagram upto 7 groups
# require: library(ggVennDiagram);library(ggplot2);library(nVennR)
#install.packages('/storage/singlecell/sangbaek/software/nVennR_0.2.3.tar.gz')
# user parameters: gList-id list per group,prefix, catNames-category name(group names)
# usage: venPlot <- runVennDgramG7(gList=x,prefix)  
runVennDgramG7 <-function(gList,prefix,catNames=NULL){
  if(is.null(catNames)){catNames <- names(gList)}else{catNames} # replace with new names}
  grnumb   <- length(catNames)
  
  # plot
  p <- ggVennDiagram(gList,category.names = catNames,label_color = "black", label_size = 4)# + scale_fill_gradient(low="blue",high = "red")
  p2<- p + scale_fill_distiller(palette = "RdBu")
  pdf(paste0("plot_venn_g",grnumb,"_",prefix,".pdf"))
  print(p2)
  dev.off()
  
  # get classified list
  #overlap <- calculate.overlap(xx.1)
  overlap <- calculate.overlap(gList)
  
  #save 
  saveRDS(overlap, paste0('obj_overlap_summary_venn_table_groups_',grnumb,'_',prefix,'.rds'))
  sink(paste0('summary_venn_table_groups_',grnumb,'_',prefix,'.txt'))
  print(str(overlap))
  print(overlap)
  sink()
  
  #library(nVennR)
  myV <- plotVenn(gList, showPlot = F, outFile = paste0("plot_venn_g",grnumb,"_",prefix,"_bynVennR.svg"))
  sink(paste0('summary_venn_table_groups_',grnumb,'_',prefix,'_bynVenn.txt'))
  print(str(myV))
  print(listVennRegions(myV))
  sink()
  
  return(p2)
}



# function: get GRange object
getGRanges<-function(dirIn, finbed){
  #setwd(dirIn)
  peaksA <- read.table(file = finbed,col.names = c("chr", "start", "end"))
  gr.A   <- makeGRangesFromDataFrame(peaksA)
  return(gr.A)
}



# Function-get ATAC signac obj
# parameters: finObj-sigobj_qc,finfrag-fragmentFile,annotations=annotObject,prefix,combined.peaks
# usage:  aAtac <- getEachAtacObjwcomFeatures(finObj,finfrag,annotations,prefix,combined.peaks)

getEachAtacObjwcomFeatures <-function(finObj,finfrag,annotations,prefix,combined.peaks){
  atac1<- readRDS(finObj)
  md.A <- atac1@meta.data
  DefaultAssay(object = atac1) <- "peaks"
  finfrag <- paste0(dirRawList[i],'/fragments.tsv.gz')
  frags.A <- CreateFragmentObject(path=finfrag, cells=colnames(atac1))
  
  atac1.counts <- FeatureMatrix(fragments = frags.A, features = combined.peaks,
                                cells = colnames(atac1))
  
  atac1_assay <- CreateChromatinAssay(atac1.counts, fragments = frags.A)
  
  aAtac1      <- CreateSeuratObject(atac1_assay, assay = "ATAC", meta.data=md.A)
  Annotation(aAtac1) <- annotations
  aAtac1$dataset <- prefix
  print(head(colnames(aAtac1)))
  # save RDS
  saveRDS(aAtac1, paste0("obj_annotAtac_",prefix,"_wcompeaks.rds"))
  return(aAtac1)
}




# function-merge atac signac objects
# parameters: aAtacList-ATAC obj list,prefixList
# merge all datasets, adding a cell ID to make sure cell names are unique
# usage: combined <- getCombinedSignacObjects(aAtacList,prefixList)
getCombinedSignacObjects<- function(aAtacList,prefixList){
  combined <- merge(
    x = aAtacList[[1]],
    y = aAtacList[2:length(aAtacList)],
    add.cell.ids = prefixList )
  #combined[["ATAC"]]
  print(combined)
  
  # clustering
  combined <- RunTFIDF(combined)
  combined <- FindTopFeatures(combined, min.cutoff = 20)
  combined <- RunSVD(combined)
  combined <- RunUMAP(combined, dims = c(2:30), reduction = 'lsi')
  combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = c(2:30))
  combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3)
  return(combined)
}



# function-Assign cell types
# parameters: objectA, prefix,dirOut,newCellTypes-new celltype list
# usage: objA <- runAssignCelltype(objectA, prefix, dirOut,newCellTypes)

### Assign cell types ------
runAssignCelltype <- function(objectA, prefix, dirOut,newCellTypes){
  names(newCellTypes) <- levels(objectA)
  objectA <- RenameIdents(objectA, newCellTypes)
  objectA$celltype <- Idents(objectA)
  
  # umap
  #  if(F){
  pdf(paste0("plot_umap_d20_wsubset_wCelltypes_",prefix,".pdf"))
  p1 <- DimPlot(object = objectA, reduction = "umap",group.by="celltype",label=T,label.size = 8)
  p2 <- DimPlot(objectA,reduction = "umap",group.by="subset",label=F,label.size = 5)
  p3 <- DimPlot(objectA,reduction = "umap",group.by="subset",split.by = "subset",label=F,label.size = 5)
  p4 <- DimPlot(objectA,reduction = "umap",group.by="celltype",split.by="subset",label=F,label.size = 5)
  print(p1)
  print(p2)
  print(p3/p4)
  dev.off()
  #  }
  return(objectA)
}

### Reclustering for Seurat object of scRNA data ----
runReclust <- function(objectA,prefix, dirOut,rsN,dn){
  # parameters to tiflter
  #qtmin    <- 0.03 # lower quantile for cut-off 
  #qtmax    <- 0.99 # higher quantile for cut-off 
  if(F){
    if(is.null(qtmin)) {qtmin <- 0.03} else {qtmin <- qtmin}
    if(is.null(qtmax)) {qtmax <- 0.99} else {qtmax <- qtmax}
    if(is.null(rsN)) {rsN <- 0.3} else {rsN <- rsN}
    if(is.null(cutMT)) {cutMT <- 20} else {cutMT <- cutMT}
    if(is.null(normOpt)) {normOpt <- "SCT"} else {normOpt <- normOpt}
    if(is.null(dn)) {dn <- 20} else {dn <- dn}
  }
  
  # set dir
  subDir <- paste0("out_reclustering_",prefix)
  setwd(dirOut)
  dir.create(subDir)
  setwd(subDir)
  
  
  # plot format: p02 +p02 + plot_layout(ncol = 2, width = c(1, 2))
  pdf(paste0("1.1.plot_QC_box_filter_bycluster_scRNA_after.filted_",prefix,".pdf"), width = 12,height = 6)
  p1 <- VlnPlot(object = objectA,pt.size = 0, features = c("nFeature_RNA","nCount_RNA","percent.mito"),group.by = "subset", ncol = 3)
  p2 <- VlnPlot(object = objectA,pt.size = 0, features = c("nFeature_RNA","nCount_RNA","percent.mito"), ncol = 3)
  print(p1/p2)
  dev.off()
  
  objectA <- RunPCA(object = objectA, features = VariableFeatures(object = objectA))
  
  #Run non-linear dimensional reduction (UMAP/tSNE)
  objectA <- RunUMAP(object = objectA, reduction = "pca", dims = 1:dn)
  
  # cell clustering
  objectA <- FindNeighbors(object = objectA, reduction = "pca", dims = 1:dn)
  objectA <- FindClusters(objectA, resolution = rsN)
  
  # umap
  pdf(paste0("1.2_plot_umap_d20_rs_",rsN,"_",prefix,".pdf"))
  print(DimPlot(object = objectA, reduction = "umap",label=T,label.size = 8))
  dev.off()
  
  cat("Process was done!\n")
  return(objectA)
  
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
  aTop2  <- objA.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
  aTop5  <- objA.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
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
# test DEGs
# runDEGs(objectA,prefix = "test",dirOut = dirOut


checkMarkersHsHasATAC <- function(objectA,markerSet,prefix,dirOut){
  #markerSet #Mlist1
  markerSetNames <- names(markerSet)
  print(names(markerSet))
  # set dir
  subDir <- paste0("out_check_markerGeneExpression_",prefix)
  setwd(dirOut)
  dir.create(subDir)
  setwd(subDir)
  
  allGenes <- rownames(objectA)
  for (i in 1:length(markerSetNames)) {
    cat ("process_marker_set: ",i,"\n")
    mgene1   <- markerSet[[i]]
    mgeneIS1 <- toupper(mgene1)
    # set plot size
    iwid <- 8+length(mgeneIS1)*1/3
    if(length(mgeneIS1)/4 <=1){iwid2<-300*length(mgeneIS1)} else{iwid2= 1200}
    ihig2=300*ceiling((length(mgeneIS1)/4))
    ihig3=3*ceiling((length(mgeneIS1)/4))  
    # dot plot
    pdf(paste0(i,"_1_plot_dot_markers_",markerSetNames[i],"_",prefix,".pdf"),width = iwid)
    print(DotPlot(object = objectA,features = unique(mgeneIS1), dot.scale = 6, assay = "ATAC",col.min = 0, scale.max = 80,
                  cols = c("gray","magenta1")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
            theme(panel.grid.major = element_line(colour = "grey95", size = (.1)))+ggtitle(paste0("Marker gene expression:",markerSetNames[i],"- ",prefix))) #+ coord_flip
    dev.off()
    
    # feature plot
    #pdf(paste0(prefix, "_2.2_featureplot_celltype.markers.pdf"),width = 9, height =15 )
    png(paste0(i,"_2_featureplot_celltype_markers_",markerSetNames[i],"_",prefix,"2.png"),width = iwid2, height = ihig2)
    print(FeaturePlot(object = objectA, features = mgeneIS1, pt.size=.1, min.cutoff = 0, ncol = 3)+
            labs(caption = paste0("TissueMarker: ",markerSetNames[i],"- Sample: ",prefix)))
    dev.off()
    # vlnplot
    #png(paste0(prefix,"_2_4_vlnplot_celltype.markers_",markerSetNames[i],"_",prefix,".png"),width = iwid2*1.5, height = ihig2))
    pdf(paste0(i,"_3_vlnplot_celltype_markers_",markerSetNames[i],"_",prefix,".pdf"),width = 12, height =ihig3)
    print(VlnPlot(object = objectA,features = mgeneIS1,pt.size = 0, ncol = 3)+
            labs(caption = paste0("TissueMarker: ",markerSetNames[i],"- Sample: ",prefix)))
    dev.off()
  }
  #print(markerSet)
  print("Process was done!")
  
  tg <-c("Rho","Arr3","Opn1mw","Glul","Grm6","Prkca","Vstm2b","Grik1","Sox6","Slitrk5","Gad1","Calb2","Thy1","Lhx1","C1qa","Rpe65")
  tg <- toupper(tg)
  #"Grm6(BC.on)","Prkca(BC)","Vstm2b(RBC)","Grik1(BC.off)","Sox6(BC.5A)","Slitrk5(BC.5C)","Gad1","Calb2(AC.sac)","Thy1(RGC)","Lhx1(HC)","C1qa(Mic)"
  ### check marker-feature/Vlnplot
  pdf(paste0(prefix,"_2_1_dotplot_celltype.marker.pdf"),width = 12, height =9 )
  print(DotPlot(object = objectA,features = tg, dot.scale = 6, assay = "RNA",col.min = 0, scale.max = 80,
                cols = c("gray","magenta1")) + RotatedAxis() + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
          theme(panel.grid.major = element_line(colour = "grey95", size = (.1)))) #+ coord_flip
  dev.off()
  
  #pdf(paste0(prefix, "_2.2_featureplot_celltype.markers.pdf"),width = 9, height =15 )
  png(paste0(prefix,"_2_2_featureplot_celltype.markers.png"), width=800, height = 1200, units="px")
  print(FeaturePlot(object = objectA, features = tg, pt.size=.1, min.cutoff = 0, ncol = 3))
  dev.off()
  
  pdf(paste0(prefix,"_2_3_featureplot_celltype.markers.pdf"), width=9, height = 3)
  #png(paste0(prefix,"_2.3_featureplot_celltype.markers.png"), width=600, height = 200, units="px")
  print(FeaturePlot(object = objectA,pt.size = 0, features = c("nFeature_RNA","nCount_RNA","percent.mito"), ncol = 3))
  dev.off()
  
  pdf(paste0(prefix, "_2_4_vlnplot_celltype.markers.pdf"),width = 12, height =15 )
  print(VlnPlot(object = objectA,features = tg,pt.size = 0, ncol = 3))
  dev.off()
}

# function - generate multiple umap with 6 different resolution (0.1 - 1)
# library: Seurat, cowplot,ggplot2,dplyr
# parameters: objectA,prefix,dm
# usage: objA <- checkUmapWmultipleResolution(objectA,prefix,dm)

checkUmapWmultipleResolution <- function(objectA,prefix,dm=20){
  # clustering
  objectA <- FindNeighbors(objectA, dims = 1:dm)
  objectA <- FindClusters(objectA, resolution = c(0.1,0.2,0.4, 0.6, 0.8,1), save.SNN = TRUE)  
  objectA <- RunUMAP(objectA, dims = 1:dm) 
  names(head(objectA))
  
  # plot
  png(paste0("umap_multi_rs_n6_",prefix,".png"), width=1400, height = 800, units="px")  
  plot_grid(nrow=2, ncol = 3, DimPlot(objectA, reduction = "umap", group.by = "RNA_snn_res.0.1") +
              ggtitle("RNA_snn_res.0.1"), DimPlot(objectA, reduction = "umap", group.by = "RNA_snn_res.0.2") +
              ggtitle("RNA_snn_res.0.2"), DimPlot(objectA, reduction = "umap", group.by = "RNA_snn_res.0.4") +
              ggtitle("RNA_snn_res.0.4"), DimPlot(objectA, reduction = "umap", group.by = "RNA_snn_res.0.6") +
              ggtitle("RNA_snn_res.0.6"), DimPlot(objectA, reduction = "umap", group.by = "RNA_snn_res.0.8") +
              ggtitle("RNA_snn_res.0.8"), DimPlot(objectA, reduction = "umap", group.by = "RNA_snn_res.1") +
              ggtitle("RNA_snn_res.1"))
  dev.off()
  
  #save RDS
  saveRDS(objectA,paste0("obj_umap_multi_rs_n6_",prefix,"_dm_",dm,".rds") )
  
  return(objectA)
}



### Functtion- draw marker plots: dot, vln and feature plots
# reauired lib: Seurat, ggplot2
# objectA, mkers- marker list, prefix
# usage: DrawMker3Plot(objectA, mkers, prefix)

DrawMker3Plot<- function(objectA, mkers, prefix){
  cat("the dir for ouput plots: ",getwd(),"\n")
  #print(mode(mkers))
  prefix <- paste0(prefix,"_knownM")
  ya<-round((length(mkers)/4))*300
  xl<-round((length(mkers)/7))+5
  y1<-round(length(table(Idents(objectA)))/6)+5
  
  #umap
  p1 <- DimPlot(objectA, label = TRUE, repel = TRUE)
  
  # dotplot
  p2<- DotPlot(object = objectA,features = mkers, dot.scale = xl, assay = "RNA",col.min = 0, 
               scale.max = 80, cols = c("gray","magenta1")) + RotatedAxis() + 
    theme(axis.text.x = element_text(angle = 45,hjust=1)) +
    theme(panel.grid.major = element_line(colour = "grey95", size = (.1))) #+ coord_flip
  
  pdf(paste0("01_plot_dot_markers","_",prefix,".pdf"), width = xl, height = y1)
    print(p2 + ggtitle(paste0("Dotplot-",prefix))) 
  dev.off()
  
  p3 <-VlnPlot(object = objectA, features = mkers, pt.size=0, ncol = 4)
  png(paste0("02_plot_vlin_",prefix,".png"),width = 1200, height = (ya+300))
    print(p3+p1+ ggtitle(paste0("UMAP-",prefix)))
  dev.off()
  
  p4 <- FeaturePlot(objectA,features = mkers,pt.size=0, ncol = 4)
  png(paste0("03","_plot_feature_",prefix,".png"),width = 1200, height =(ya+300))
    print(p4+ p1+ ggtitle(paste0("UMAP",prefix)))
  dev.off()
}


# Function: get HAS marker List for "CiliaryBody","Cornea","Iris","Lens","Limbus", "TM"
# usage: mkerHas <- getHasMarkers()
getHasMarkers <- function() {
  markerSet <- readRDS("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/marker_genes_HAS/obj_marker_HAS_list_6T.rds")
  print(markerSet)
  return(markerSet)
}


## function - convert rmd to r script
# parameter: input-rmd file
# required R pkg: rmarkdown, knitr package
# usage: convert_ipynb2RmdRscript(input)
convert_ipynb2RmdRscript <-function(input){
  output1 <- gsub(".ipynb",".Rmd",input)
  output2 <- gsub(".ipynb",".R",input)
  rmarkdown ::convert_ipynb(input, output1)
  knitr::purl(output1, output2, documentation = 2)
}



