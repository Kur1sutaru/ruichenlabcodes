library(Seurat)
library(SeuratDisk)
#library(rhdf5)
args <- commandArgs(trailingOnly = TRUE)
sam=args[1]
dir="/storage/chentemp/wangj/human_ret_anc/data/HRCA/"
Convert(paste0(dir ,sam,".h5ad"), dest = "h5seurat", overwrite = TRUE)

seRNA0 <- LoadH5Seurat(paste0(dir,sam,".h5seurat"),misc = FALSE,meta.data = FALSE)

meta.data=read.table(paste0(dir, sam,".obs.txt.gz"),header=T,sep="\t")

#rownames(meta.data)=meta.data$barcode

seRNA=AddMetaData(object=seRNA0, metadata=meta.data)

saveRDS(seRNA,file=paste0(dir, sam,".rds"))

