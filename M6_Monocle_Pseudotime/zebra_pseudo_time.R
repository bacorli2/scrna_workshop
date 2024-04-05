

library(tidyverse)
library(Seurat)
library(patchwork)
library(HGNChelper)
library(openxlsx)
library(presto)
library(scAnnotatR)
library(SingleR)
library(celldex)
library(monocle3)
library(SeuratWrappers)
library(cowplot)
require(DESeq2)
# library(DESeq2)
options(ggrepel.max.overlaps = Inf) 

# Set wd to base of workshop repository
here::i_am("README.md")
source(here::here("R_override", "scType_SeuratObj.R"))

# 10X Dataset Importation ################################################
#_______________________________________________________________________________

pca_dims <- 1:20

data_paths = c("Day_5", "Day_10")
srat_list = list()
for (n in 1:2) {
  # Load count, cell, and gene data files
  data_10x <- Read10X(data.dir = here::here("_temp_data", "H5AD", data_paths[[n]]))
  srat_list[[n]] <- CreateSeuratObject(counts = data_10x, project = "zebra", 
                                      min.cells = 3, min.features = 200)
  # Add additional cell type labels
  meta.data <- read.csv(here::here("_temp_data", "H5AD", data_paths[[n]], "metadata.csv"))
  srat_list[[n]] <- AddMetaData(srat_list[[n]], meta.data)
  
}


# Seurat Merge, Initial Processing #############################################
#_______________________________________________________________________________

# Merge seurat datasets into one object (does not integrate across datasets)
srat.comb <- merge(srat_list[[1]], y = srat_list[[2]], add.cell.ids = c("5", "10"), 
                   project = "Zebra")
# Initially process data without integration between datasets.
srat.comb <- NormalizeData(srat.comb)
srat.comb <- FindVariableFeatures(srat.comb)
srat.comb <- ScaleData(srat.comb)
srat.comb <- RunPCA(srat.comb)
ElbowPlot(srat.comb, ndims = 30)
srat.comb <- FindNeighbors(srat.comb, dims = pca_dims, reduction = "pca")
srat.comb <- FindClusters(srat.comb, resolution = 2, 
                               cluster.name = "unintegrated_clusters")
srat.comb <- RunUMAP(srat.comb, dims = pca_dims, reduction = "pca", 
                          reduction.name = "umap.unintegrated")


# Integrate Datasets, Reprocess ################################################
#_______________________________________________________________________________

srat.comb <- IntegrateLayers(object = srat.comb, method = CCAIntegration, 
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# Merge datasets into a single layer now that we have integrated
srat.comb <- JoinLayers(srat.comb)
srat.comb <- FindNeighbors(srat.comb, reduction = "integrated.cca", dims = pca_dims)
srat.comb <- FindClusters(srat.comb, resolution = .5, cluster.name = "cca_clusters")
srat.comb <- RunUMAP(srat.comb, reduction = "integrated.cca", dims = pca_dims, 
               reduction.name = "UMAP")
# DimPlot(srat.comb, reduction = "umap.cca", group.by = c("Method"),
#   combine = FALSE, label.size = 2)
srat.comb$seurat_clusters = srat.comb$cca_clusters

DimPlot(srat.comb, reduction = "UMAP", group.by = "cca_clusters",
        combine = FALSE, label.size = 2)

DimPlot(srat.comb, reduction = "UMAP", group.by = "seurat_clusters",
        combine = FALSE, label.size = 2)



# Convert to Monocle ##########################################################
#_______________________________________________________________________________

# Convert seurat object to monocle with wrapper (srat object from beginning)
cds <- SeuratWrappers::as.cell_data_set(srat.comb)

## Calculate size factors using built-in function in monocle3, add gene names
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- 
  rownames(srat.comb[["RNA"]])


# Reclustering #################################################################
#_______________________________________________________________________________

# Monocle needs partitions as well as clusters
# Using cluster_method = leiden raises error with Nonsymmetric adjacency matrix
# cluster_method = c("leiden", "louvain")
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 20, num_iter = 1,
                     cluster_method = "leiden", partition_qval = 0.05,
                     resolution = 1.5e-3,
                     weight = FALSE, random_seed = 1,  verbose = FALSE)



# 3A) Monocle needs cells groups by partition (metaclusters) as well as clusters
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE, 
           group_label_size = 4)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE,
           group_label_size = 4)
# plot_cells(cds, color_cells_by = "cell_type", show_trajectory_graph = FALSE,
#            group_label_size = 4, alpha = 0.5)



# Find Cell Trajectories #######################################################
#_______________________________________________________________________________

# 4) Fit principal graph (trajectory) for each partition
# Fits trajectories between projected cells
# Identifies branchpoints and endpoints with cell trajectories
cds <- learn_graph(cds)

plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE,
           label_leaves = TRUE, label_principal_points = TRUE, 
           graph_label_size = 3, alpha = 0.5)


# NOTE: The node ID will chang between systems, you need to look at the graph 
# above and change the node id "Y_xxx" and "Y_xxx" in the function calls below.

# 5) Order cells in pseudotime from selected principle node(s)
cds <- order_cells(cds, root_pr_nodes = c("Y_91"))
# 5A) visualize pseudotime from root node
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, 
           graph_label_size = 1.5, alpha = 0.5)


# Plot this and choose points for the next part ################################
#_______________________________________________________________________________
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE,
           label_leaves = TRUE, label_principal_points = TRUE, 
           graph_label_size = 3, alpha = 0.5)
# Find an interesting branch node and end node you want to subset, 
# and replace these two arguments in the next step.
# starting_pr_node
# ending_pr_nodes


# Choose Specific Trajectory ###################################################
#_______________________________________________________________________________

# 6) Subset cells from a particular trajectory (graph segment)
cds_sub <- choose_graph_segments(cds, reduction_method = "UMAP",
                                 starting_pr_node = "Y_231", 
                                 ending_pr_nodes = c("Y_281"),
                                 clear_cds = FALSE) 


# Keep previous normalization, so we proceed straight to reclustering
# # 3) Cluster cells in reduced space
cds_sub <- cluster_cells(cds_sub,reduction_method = "UMAP", k=100)
# 4) Construct new graph
cds_sub <- learn_graph(cds_sub)
# 4A) Color cells by celltype and and label principle points
plot_cells(cds_sub, color_cells_by = "cluster", 
           label_cell_groups = FALSE, label_leaves = TRUE, 
           label_principal_points = TRUE, alpha = 0.5,
           graph_label_size = 3)

# 5) Genes that are associated with the selected trajectory 
#  (cores = 1) for reproducibility
subset_pr_test_res <- graph_test(cds_sub, neighbor_graph="principal_graph", 
                                 cores = 1)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))


# Find Gene Modules ############################################################
#_______________________________________________________________________________

# 6) Group genes into modules to visualize expression trends over pseudotime
cds_sub <- preprocess_cds(cds_sub, num_dim = 15)
gene_module_df <- find_gene_modules(cds_sub[pr_deg_ids,])

# 7) Order modules by similarity (via hclust) to see which ones activate earlier
agg_mat <- aggregate_gene_expression(cds_sub, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)
                                [module_dendro$order])

# Visualize gene module activation over pseudotime
# rowData(cds_sub)$gene_short_name <- rowData(cds)$gene_name
plot_cells(cds_sub, genes=gene_module_df)




