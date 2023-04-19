setwd("D:/baylor/TRIGEMINAL/B6_TG_Nu/analysis/merfishpanelgenes")
library(Seurat)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)
library(mclust)
# load data

objListallgenes <- Read10X_h5(filename = "D:/baylor/TRIGEMINAL/B6_TG_Nu/analysis/input/filtered_feature_bc_matrix.h5")

objListallgenes <- CreateSeuratObject(counts = objListallgenes, project = "B6_TG_Nu", min.cells = 3, min.features = 200)
objList

dense.size <- object.size(objListallgenes)
dense.size

sparse.size <- object.size(objListallgenes)
sparse.size

objListallgenes[["percent.mt"]] <- PercentageFeatureSet(objListallgenes, pattern = "^MT-|^Mt-")

objListallgenes <- NormalizeData(objListallgenes, normalization.method = "LogNormalize", scale.factor = 10000)

objListallgenes <- FindVariableFeatures(objListallgenes, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(objListallgenes), 10)
write.csv(top10,"B6_TG_Nutophighvariablegenes.csv")

all.genes <- rownames(objListallgenes)
objListmerfish <- ScaleData(objListallgenes, features = c("Cngb1","Pde6a","Reep6","Nt5e","Susd3","Kcne2","Adrb1","Cngb3","Pde6c","Casp7","Id3","Kdr","Vim","Rdh10","Kcnj10","Onecut2","Onecut1","Lhx1","Lima1","Adra2a","Otx2","Vsx2","Grm6","Isl1","Grik1","Crx","Vsx1","Pcdh17","Mylk","Tacr3","Pcdh10","Nxph1","Wls","Ncam2","Ebf1","Pde1a","Syt2","Chrnb3","Erbb4","Irx5","Chrna6","Irx6","Gng2","Slc16a7","Lrrtm3","Cacna2d2","Col11a1","Calml4","Nxph4","Sox6","Ptprt","Pla2g7","Chrm2","Ddit4l","Gas1","Nwd2","Areg","Slitrk5","Msrb2","Pcdh8","Lrrtm1","Cplx2","Irx3","Kirrel3","Spock3","Cdc25c","Reln","Igfn1","Rgs3","Kcnab1","Sox5","Cpne9","Serpini1","Plcb1","Spon1","Prkca","Strip2","Tpbg","Vstm2b","Pax6","Slc32a1","Tfap2a","Gad2","Slc6a9","Slc6a1","Ptn","Bhlhe22","Crybb3","Maf","Mmp9","Gjd2","Drd2","Prom1","Igf1","Cbfa2t3","Pantr1","Crot","Nefh","Olfm3","Tmem114","Arap2","Tfap2e","Ntrk1","Trhde","Shisa9","Cntn6","Nckap5","Slc35d3","Syt10","Fgf18","Slc17a8","Gulo","Igfbp5","Vmn2r122","Wnt5a","Rnf152","Gng7","Slc18a3","Cpne4","Mfap5","Nxph2","Gdf3","Sema3a","Rxfp3","Elavl2","C1ql3","Cbln4","Plscr1","Rasgrp1","Pdgfra","Neurod6","Col27a1","Crym","Abi3bp","Igsf3","Rab3b","Nfib","Ctsh","Glra1","Lgr5","Klc3","Chrdl1","Fgf7","Ramp3","Rspo4","Robo3","Frzb","Tpm2","Bmp3","Sgcz","Kctd12","Hmcn1","Crh","Pthlh","Bcl11b","Kcnd2","Mef2c","Akap12","Gch1","Drd5","Pou3f1","Etv1","Mlip","Tbx2","Pcdh18","Fst","Otor","Chl1","Nefm","Slc4a4","Slc38a5","Fezf2","Adra1b","Eomes","Sh3bp5","Postn","Ldb2","Gm12371","Bmp2","Ngf","Pcsk9","Tmem163","Arpp21","Mafb","Spp1","Lrtm1","Slc17a7","Tmeff2","Apela","Antxr2","Penk","Pdyn","Tnfsf8","Epha3","Cdh20","Pappa2","Cpm","Crip2","Slc17a6","Rbpms","Pou4f2","Serpine2","Amigo2","Lypd1","Foxp2","Irx4","Gpr83","Tbr1","Pcdh20","Zic1","Tbx20","Tagln2","Prkcq","Tac1","Slc7a11","Plpp4","Whrn","Serpinb1b","Mmp17","Trhr","Runx1","Vit","Meis2","Tfap2d","Col25a1","Gabrg1","4833423E24Rik","Lmo2","Prdm8","Kctd4","Prlr","Anxa3","Calca","Gpc5","Cdhr1","Prokr1","Ttn","Cxcl13","Necab1","Prkcg","Sfrp2","Gpr101","Cdk15","Stxbp6","Ppp1r1c","Ifitm10","Gas7","Kcna5","Adcyap1","Opn4","Hapln1","Igfbp4","Cbln2","Coch","Ceacam10","Rarres1","Amhr2","Neurod2","S100b","Cd24a","Ngfr","Kit","Scn4b","Fes","Il1rapl2","Tnni3","Fxyd6","Kcnip2","Hpse","Tulp1","Neurod1","Thrb","Lyar","Nrl","Nr2e3","Ybx3","Casz1","Tax1bp1","Rora","Jun","Hes1","Fos","Egr1","Junb","Zfp36l1","Tsc22d4","Id2","Sox2","Hopx","Sox9","Fosb","Lhx2","Rax","Nr2e1","Hes5","Tox3","Sox8","Jund","Creb3l2","Tox2","Tfap2b","Nhlh2","Myt1l","Zfhx3","Nr4a2","Zeb2","Tcf4","Nfia","Six6","Tshz2","Lhx9","Pbx1","Pou6f2","Prox1","Bcl11a","Neurod4","Peg3","Tshz1","Esrrg","Id4","Nr2f1","Nr2f2","Sox4","Tbx3","Zfp503","Lhx4","Zfhx4","Dach1","Hlf","Lhx3","St18","Zbtb38","Camta1","Nfix","Ebf3","Pou4f3","Satb1","Nfe2l1","Scrt1","Isl2","Pura","Dmrtb1","Foxp1","Purb","Barhl2","Bhlhe40","Dbp","Klf13","Onecut3","Tsc22d3","Aqp6","Drd4","Nell1","Endou","Cpne7","Ptx3","Nos1","Kcnk15","Lmo3","Cacna1e","Igfbp2","Gldn","Calb2","Slc1a3","Cabp5","Hcn1","Tppp3","Vip","Atp2b1","Sgk1","Rcvrn","Hif1a","Scgn","Cacna2d3","Cnmd","Clu","Gabrr2","Gnb3","Prph2","Syt4","Pdcd4","Pde6g","Rp1","Atp1b1","Nebl","Chgb","Mcub","Calb1","Gnb1","Gngt1","Rho","Dkk3","Tent5a","Crabp1","Rbp3","C1ql2","Slc24a1","Patj","Car10","Pde3a","Pde1c","Cplx3","Tgfb2","Zfp804a","Scg2","Lrp1b","Ppp1r17","Ndnf","Rbfox1","Lamp5","C1ql1","Cntnap2","Pcp4","Trarg1","Adarb2","Zfp365","Gm4792","Car8","Gabrg3","Lrrtm4","Kcnh7","Sag","Eno1","Csmd1","Glul","Prune2","Gsg1","Pde10a","Tafa4","Grm8","Gpr179","Necab2","Slc12a5","Rprm","Opn1mw","Tenm2","Grik2","Camk2d","Dab1","Slc24a2","Ptprd","Gngt2","Gria2","Nnat","1-Mar","Pcp4l1","Rs1","Rorb","Cntn4","Arr3","Dmd","Tmsb10","Cadm1","Zbtb20","Ccdc136","Lrfn5","Egfem1","Lmo4","Grm5","Trnp1","Cadm2","Ntng1","Pde6h","Slc24a3","Pou4f1","Gad1","Cdh8","Pcp2","Acsl3","Kcnma1","Lsamp","Opn1sw","Cntn5","Csmd3","Pdc","Apoe","Unc5d","Gng13","Slit2","Kcnb2","Six3","Kcnip4","Rlbp1","Rpgrip1","Ank3","Unc13c","Kcnq1ot1","Rgs6","Rd3l","Nlgn1","Pcdh7","Pcdh15","Ncam1","Pcdh9","Nrxn1","Zfp804b","Gabra1","Trpm1","Ptprk","Trpm3"))
objListmerfish <- RunPCA(objListmerfish, features = VariableFeatures(object = objListmerfish))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(objListmerfish), 10)
write.csv(top10,"B6_TG_NuobjListmerfishtophighvariablegenes.csv")

# plot variable features with and without labels
jpeg(file="B6_TG_NuobjListmerfishvariable_plot1.jpeg")
plot1 <- VariableFeaturePlot(objListmerfish)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

objListmerfish <- FindClusters(objListmerfish, resolution = 0.5)
jpeg(file="objListmerfishumapres0.5.jpeg")
DimPlot(objListmerfish, reduction = "umap")
dev.off()

jpeg(file="B6_TG_Numerfish.jpeg")
DimPlot(objListmerfish, reduction = "pca")
dev.off()

jpeg(file="B6_TG_NuobjListmerfishpcascorepc1.jpeg")
DimHeatmap(objListmerfish, dims = 1, cells = 500, balanced = TRUE)
dev.off()

jpeg(file="B6_TG_NuobjListmerfishpcasdim15.jpeg")
DimHeatmap(objListmerfish, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

objListmerfish <- JackStraw(objListmerfish, num.replicate = 100)
objListmerfish <- ScoreJackStraw(objListmerfish, dims = 1:20)

jpeg(file="B6_TG_NuobjListmerfishJackStrawPlot30pc.jpeg")
JackStrawPlot(objListmerfish, dims = 1:20)
dev.off()

jpeg(file="B6_TG_NuobjListmerfishElbowPlot.jpeg")
ElbowPlot(objListmerfish)
dev.off()


objListmerfish <- FindNeighbors(objListmerfish, dims = 1:10)
objListmerfish <- FindClusters(objListmerfish, resolution = 0.5)

objListmerfish <- RunUMAP(objListmerfish, dims = 1:10)
jpeg(file="objmerfish0.510dims.jpeg")
DimPlot(objListmerfish, reduction = "umap", label=T)
dev.off()

# Look at cluster IDs of the first 5 cells
identsmerfish <- head(Idents(objListmerfish), 5)
write.csv(identsmerfish, "identsmerfish.csv")

identsmerfish <- head(objListmerfish$orig.ident, 5)
write.csv(identsmerfish, "identsidentsmerfish.csv")


jpeg(file="objListmerfishobjListmerfishumapres0.5.jpeg")
DimPlot(objListmerfish, reduction = "umap")
dev.off()

objListmerfish <- FindClusters(objListmerfish, resolution = 0.3)
jpeg(file="objListmerfishumapres0.3.jpeg")
DimPlot(objListmerfish, reduction = "umap")
dev.off()

objListmerfish <- FindClusters(objListmerfish, resolution = 1)
jpeg(file="objListmerfishumapres1.jpeg")
DimPlot(objListmerfish, reduction = "umap")
dev.off()

objList <- FindClusters(objList, resolution = 0.4)
jpeg(file="objListmerfishumapres0.4.jpeg")
DimPlot(objListmerfish, reduction = "umap")
dev.off()

# find all markers of cluster 2
cluster2.markers <- FindMarkers(objListmerfish, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
write.csv(cluster2.markers, "merfishgenescluster2.markers.csv")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
objList.markers <- FindAllMarkers(objListmerfish, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
objList.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(objList.markers, "objListmerfishobjListclustermarkers.csv")


objList.markers <- FindAllMarkers(objListmerfish, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
objList.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(objList.markers, "objList.markersmerfish.csv")

cluster0.markers <- FindMarkers(objListmerfish, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


jpeg(file="objListmerfishobjListtestviolinfeatplot.jpeg")
VlnPlot(objListmerfish, features = c("Sulf1"))
dev.off()

jpeg(file="objListmerfishobjListtesfeatplot0.jpeg")
FeaturePlot(objListmerfish, features = c("Sulf1","Aff3","Pid1","Sned1","Susd3","Steap3"))
dev.off()

jpeg(file="objListmarkersheatmap.jpeg")
objList.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(objList, features = objListmerfish) + NoLegend()
dev.off()
