##Retrieve the bulk counts from single cell
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(hdf5r)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DESeq2)
library(rtracklayer)

bulk_data <- AverageExpression(onhannotated)
write.csv(bulk_data$RNA, "onhbulk.csv")

#Tranform into FPKM values
bulk_counts <- read.csv("bulk_expression.csv", row.names = 1)
dds <- DESeqDataSetFromMatrix(countData = bulk_counts, colData = colData, design = ~ condition)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

#### Calculate Gene Length ###
# You'll need information about the length of each gene. 
# If you don't have this information, you may need to obtain 
# a reference genome annotation file (GTF/GFF) and extract the gene lengths.

#How to get a GTF file
# Create an AnnotationHub object
ah <- AnnotationHub()

# Search for mouse GTF files
mouse_gtf <- query(ah, c("mus musculus", "gtf"))

# Print the available GTF files for mouse
mouse_gtf


# Example if you have a GTF file

mouse_gtf <- import.gff("path/to/your/annotation.gtf")
gene_lengths <- width(mouse_gtf)

##Calculate FPKM values

gene_lengths_kb <- gene_lengths / 1000
fpkm_values <- normalized_counts / (gene_lengths_kb * colSums(normalized_counts) / 1e6)
