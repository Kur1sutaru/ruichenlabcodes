#To write a table with the significative p-values

setwd("D:/baylor/Dry eye/DEGS PER CLUSTER")

#To human data, install the library org.Hs.eg.db first!
#TO LOAD TH PACKAGE
library(pathfindR)



#Download and load the pathfindr package
#install.packages("pathfindR")
library(pathfindR)

pathwaysresults<-run_pathfindR(c2_resSig_clusters_cutoff0.64, 
                               gene_sets = "KEGG",
                               min_gset_size = 2,
                               max_gset_size = 300,
                               custom_genes = NULL,
                               custom_descriptions = NULL,
                               pin_name_path = "Biogrid",
                               p_val_threshold = 0.05,
                               visualize_enriched_terms = FALSE,
                               max_to_plot = NULL,
                               convert2alias = FALSE,
                               enrichment_threshold = 0.05,
                               adj_method = "bonferroni",
                               search_method = "GR",
                               plot_enrichment_chart = TRUE,
                               output_dir = "c2_resSig_clusters_cutoff0.64_pathfindR_Results",
                               list_active_snw_genes = TRUE,
                               silent_option = TRUE)

cluster_enriched_terms(pathwaysresults, use_description = TRUE)

input_results<-input_processing(
  c11_resSig_clusters_cutoff0.64,
  p_val_threshold = 0.05,
  pin_name_path = "Biogrid",
  convert2alias = TRUE
)

visualize_terms(pathwaysresults,
                input_processed = input_results,
                hsa_KEGG = FALSE,
                pin_name_path = "Biogrid")

# The table must be contain 2 or 3 columns (gene name, p value, and or fold change).
write.csv(pathwaysresults, "c11_resSig_clusters_cutoff0.64kegg.csv")
enrichment_chart(pathwaysresults, top_terms = 20)
term_gene_graph(pathwaysresults, num_terms = 3,use_description = TRUE)
term_gene_heatmap(pathwaysresults, out_top100markers_eachCluster_ciliary_combined_degs, num_terms = 5, use_description = TRUE)
