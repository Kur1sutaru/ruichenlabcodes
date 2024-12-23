##### analysis of bulk rna seq ###
## install packages if needed ##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages(c("ggplot2", "pheatmap", "EnhancedVolcano"))

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Step 1: Import count data and metadata
count_data <- read.csv("path/to/your/count_matrix.csv", row.names = 1)
metadata <- read.csv("path/to/your/metadata.csv")  # Ensure metadata has a column 'condition'

# Step 2: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = metadata,
  design = ~ condition
)

# Step 3: Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Step 4: Perform differential expression analysis
dds <- DESeq(dds)

# Step 5: Extract results
res <- results(dds, alpha = 0.05)
res <- res[order(res$padj), ]  # Sort by adjusted p-value
write.csv(as.data.frame(res), "DE_results.csv")

# Step 6: Data normalization for visualization
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "Normalized_counts.csv")

# Step 7: Plot MA plot
plotMA(res, main = "MA Plot", ylim = c(-5, 5))

# Step 8: Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "pvalue",
                title = "Volcano Plot",
                pCutoff = 0.05,
                FCcutoff = 1
)

# Step 9: Heatmap of top differentially expressed genes
top_genes <- head(order(res$padj), 20)  # Top 20 DE genes
pheatmap(
  normalized_counts[top_genes, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,
  annotation_col = metadata
)

# Step 10: PCA plot
rld <- rlog(dds, blind = TRUE)
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot")

# Step 11: Export PCA data
write.csv(pca_data, "PCA_data.csv")

# Step 12: Additional GO/Pathway analysis (optional)
# You can use clusterProfiler or similar packages for enrichment analysis if needed.

cat("Analysis completed successfully. Check the output files for results and visualizations.\n")



### Gene set enrichment analysis - GSEA
# install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))
install.packages(c("ggplot2", "tidyverse"))


# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Use org.Mm.eg.db for mouse or other databases as needed
library(enrichplot)
library(ggplot2)
library(DOSE)

# Step 1: Prepare input data
# Assuming you have a DESeq2 results table with gene names, log2FoldChange, and p-values
res <- read.csv("DE_results.csv", row.names = 1)

# Ensure gene symbols are mapped to ENTREZ IDs
gene_list <- res$log2FoldChange
names(gene_list) <- res$gene_id  # Replace 'gene_id' with the column of gene names
gene_list <- na.omit(gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

# Map gene symbols to ENTREZ IDs
gene_list_entrez <- bitr(names(gene_list), fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
gene_list <- gene_list[gene_list_entrez$SYMBOL]
names(gene_list) <- gene_list_entrez$ENTREZID

# Step 2: Perform GSEA
gsea_results <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Biological Process (can also use "MF" or "CC")
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Save GSEA results
write.csv(as.data.frame(gsea_results), "GSEA_BP_results.csv")

# Visualize GSEA results
dotplot(gsea_results, showCategory = 20) + ggtitle("GSEA Dotplot")
ridgeplot(gsea_results) + ggtitle("GSEA Ridge Plot")
emapplot(pairwise_termsim(gsea_results)) + ggtitle("GSEA Enrichment Map")

# Step 3: Perform KEGG Pathway Analysis
kegg_results <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",  # Use "mmu" for mouse, "rno" for rat, etc.
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Save KEGG results
write.csv(as.data.frame(kegg_results), "GSEA_KEGG_results.csv")

# Visualize KEGG results
dotplot(kegg_results, showCategory = 20) + ggtitle("KEGG Dotplot")
ridgeplot(kegg_results) + ggtitle("KEGG Ridge Plot")
emapplot(pairwise_termsim(kegg_results)) + ggtitle("KEGG Enrichment Map")

# Step 4: KEGG Pathway Analysis without GSEA (Optional)
# If you prefer over-representation analysis (ORA) instead of GSEA:
de_genes <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
de_genes_entrez <- bitr(de_genes$gene_id, fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)
kegg_ora <- enrichKEGG(
  gene = de_genes_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Save KEGG ORA results
write.csv(as.data.frame(kegg_ora), "KEGG_ORA_results.csv")

# Visualize KEGG ORA results
dotplot(kegg_ora, showCategory = 20) + ggtitle("KEGG ORA Dotplot")
barplot(kegg_ora, showCategory = 20) + ggtitle("KEGG ORA Bar Plot")
