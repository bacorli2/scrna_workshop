

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


# scRNA-seq Dataset Importation ################################################
#_______________________________________________________________________________
# Example small dataset (real data)
# Used from this tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial
# 2,700 single cells that were sequenced on the Illumina NextSeq 500
# 13,714 genes
# Download dataset into temp_data, unzip
pbmc3k_path <- here::here("_temp_data", "pbmc3k_filtered_gene_bc_matrices.tar.gz")
if (!file.exists(pbmc3k_path)) {
  dir.create(here::here("_temp_data"))
  download.file(paste0("https://cf.10xgenomics.com/samples/cell/pbmc3k/",
                       "pbmc3k_filtered_gene_bc_matrices.tar.gz"),
                destfile = here::here("_temp_data", 
                                      "pbmc3k_filtered_gene_bc_matrices.tar.gz"))
  untar(pbmc3k_path, exdir = here::here("_temp_data"))
}
# Load the srat dataset
srat.data <- Read10X(data.dir = here::here("_temp_data", 
                                           "filtered_gene_bc_matrices/hg19"))

# Initialize the Seurat object with the raw count matrix (non-normalized data).
# Include genes that are found within at least 3 cells
# Include cells that have at least 200 genes
srat <- CreateSeuratObject(counts = srat.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)

# Add column in metadata slot for percent of mitochondrial genes (QC metric)
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
# nFeature_RNA: total number of genes detected in each cell
# nCount_RNA: total number of molecules detected within a cell (library size)
# percent.mt: fraction of genes that are mitochondrial (qc metric)
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

## Filter poor quality cells ###################################################
# nFeature_RNA > 200: removes empty droplets or cells with little RNA
# nFeature_RNA < 25000: remove doublets (droplets with 2+ cells)
# percent.mt < 5: removes cells with over 5% mitochondrial DNA 
# (poor viability)
srat <- subset(srat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                 percent.mt < 5)


## Normalize data ##############################################################
# 1. Normalizes gene expression by the total expression in each cell
# 2. Multiplies this by a scale factor (10,000 by default)
# 3. Log-transforms the result.
# Stored in: srat[["RNA"]]$data
srat <- NormalizeData(srat, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

## Feature Selection ###########################################################
# Identify highly variables genes, (to be used for dimension reduction)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)

## Scale the data (across cells, only selected variable features) ##############
# Essentially converts gene expression to z-score (normalize by mean and std 
# across cells)
# Stored in: srat[["RNA"]]$scale.data
srat <- ScaleData(srat, features = rownames(srat))

## Scale data can also be used to remove unwanted cell cycle variation #########
# However, this is a more advance method, and it is recommended to use the new
# Seurat workflow: SCTransform(). 
# Paper: https://genomebiology.biomedcentral.com/articles/10.1186/
#        s13059-021-02584-9
# Vignette: https://satijalab.org/seurat/articles/sctransform_vignette
# srat <- ScaleData(srat, vars.to.regress = "percent.mt")

 
## Linear dimension reduction (PCA) ############################################
#_______________________________________________________________________________
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
# Plot commands: VizDimReduction(), DimPlot(), and DimHeatmap()
VizDimLoadings(srat, dims = 1:2, reduction = "pca")     
DimHeatmap(srat, dims = 1, cells = 500, balanced = TRUE)
DimPlot(srat, reduction = "pca") + NoLegend()

# Choose dimensionality of the dataset
# Maximize the signal (biological variability) to the noise (other sources of 
# variation)
ElbowPlot(srat)


## Clustering ##################################################################
#_______________________________________________________________________________
# Construct a kNN graph based on euclidean distance in a subset of PCA space 
#  (up to dimensionality chosen).
# Refine edge weights between pairs of cells based on their shared overlap and 
# local neighboors (Jaccard similarity)
srat <- FindNeighbors(srat, dims = 1:10)

# Clustering Cells: we next apply modularity optimization 
# (Louvain algorithm, SLM ) to iteratively group cells together, with the goal 
# of optimizing the standard modularity function. 
# Cluster granularity is set with resolution, 0.4-1.2 typically returns good 
# results for single-cell datasets of around 3K cells. Resolution often 
# increases for larger datasets.
srat <- FindClusters(srat, resolution = 0.5)
DimPlot(srat, reduction = "pca") + NoLegend()


# Perform UMAP clustering
srat <- RunUMAP(srat, dims= 1:10)

# Visualize UMAP clusters
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE)

# Gene expression in each cluster
# VlnPlot(data, features = c("Pax6", "Rbfox1"), slot = "counts", log = TRUE)

# FeaturePlot(data, features = c("Pax6",  "Eomes", "Aldh1l1",
#                                "Tbr1",  "Olig2", "Sox2", "Cux2", "Neurog2"))



# Dataset Integration (Simulated) ##############################################
#_______________________________________________________________________________
# For illustrative purposes, let's simulate having data from two conditions
# We can combine the data between them and instruct seurat to normalize the data
# To make comparable.
#  Groups: 0: control, 1: treatment
# We can do this by adding a factor to the seurat metadata
set.seed(0)
srat_int <- srat
srat_int@meta.data$group_id = factor(rbinom(n = ncol(srat_int), size = 1, 
                                            prob = 0.5 ), labels = c("Ctrl","Tx"))
# Split dataset based on factor column in metadata
srat_int[["RNA"]] <- split(srat_int[["RNA"]], f = srat_int$group_id)
# Integrate datasets together in seurat object
srat_int <- 
  IntegrateLayers(srat_int, method = CCAIntegration, orig.reduction = "pca", 
                  new.reduction = "integrated.cca", verbose = FALSE)
# Re-join layers after integration
srat_int[["RNA"]] <- JoinLayers(srat_int[["RNA"]])
# Rerun pipeline
srat_int <- FindNeighbors(srat_int, dims = 1:10)
srat_int <- FindClusters(srat_int, resolution = 0.5)
srat_int <- RunUMAP(srat_int, dims= 1:10)
# Visualize UMAP clusters
DimPlot(srat_int, reduction = "umap", label = TRUE,
        repel = TRUE)
# Overwrite seurat object for downstream steps
srat <- srat_int




# Cluster Marker Identification ################################################
#_______________________________________________________________________________

# Find differentially expressed genes in each cluster vs. all other clusters
# Test used is non-parametric Wilcoxon rank sum test
# Note: Install presto package for much faster results
srat_int.all.markers <- FindAllMarkers(srat_int, only.pos = TRUE)


## Cell Type Annotation: Scitype ###############################################
#_______________________________________________________________________________
#https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md
# Load gene set and cell type annotation functions into memory
source(paste0("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/",
              "master/R/gene_sets_prepare.R"))
source(paste0("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/",
              "master/R/sctype_score_.R"))
# DB file
db_ = paste0("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/",
             "master/ScTypeDB_full.xlsx")
tissue = "Immune system"
# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,
# Intestine,Muscle,Placenta,Spleen,Stomach,Thymus


# Prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# Get score matrix: cell-type (row) by cell (col)
# NOTE: scRNAseqData argument should correspond to your input scRNA-seq matrix.
#   In case Seurat is used, it is either
#   1. srat[["RNA"]]@scale.data (default),
#   2. srat[["SCT"]]@scale.data, if sctransform is used for normalization,
#   3. srat[["integrated"]]@scale.data, for joint analysis of multiple datasets.
es.max = sctype_score(scRNAseqData = srat[["RNA"]]$scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


# Merge by cluster
# For each cluster, grab all cells that below to it, find top10 best matches
# for cell type
cL_resutls = do.call("rbind", lapply(unique(srat@meta.data$seurat_clusters),
                                     function(cl){
  es.max.cl = sort(rowSums( es.max[ ,rownames(srat@meta.data[
    srat@meta.data$seurat_clusters==cl, ])]), decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl,
                  ncells = sum(srat@meta.data$seurat_clusters==cl)), 10)
}))
# Grab best cell-type match for each cluster, assign as final cell-type
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# Set low-confident (low ScType score) clusters to "Unknown"
# Sctype scores scale by n, so threshold is ncells/4
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) <
                     sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# Add column in seurat metadata for celltype annotation
srat@meta.data$cell_type <- factor(select(srat@meta.data, "seurat_clusters") %>%
  left_join(y = select(sctype_scores, "cluster", "type"),
            by = join_by(seurat_clusters == cluster)) %>% pull("type"))
# Relabel cell identity label to cell_type (previously was cluster number)
Idents(srat) <- srat@meta.data$cell_type

# UMAP Plot of Scitype annotated cells
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE,
        group.by = 'cell_type') +
  ggtitle("SciType Annotated Cells")


## Cell Classification Using scAnnotateR #######################################
#_______________________________________________________________________________
# DEFAULT MODEL: Load classification Models
default_models <- scAnnotatR::load_models("default")

# Perform classification
srat_scannot <- classify_cells(classify_obj = srat,
                             assay = 'RNA', slot = 'counts',
                             cell_types = "all",
                             path_to_models = 'default')
# Plot best match for each cell_type
DimPlot(srat_scannot, group.by = "most_probable_cell_type")



## Classification Based Cell Type: scPred ######################################
#_______________________________________________________________________________
# There is an error with scPredict and the github has not been updated in a 
# while, so we load a corrected version of function into memory.  
source(here::here("R_override", "scPredict_edited.R"))

# Process reference through seurat and scPred
ref_data <- scPred::pbmc_1 %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
ref_model <- scPred::getFeatureSpace(ref_data, "cell_type")
ref_model <- scPred::trainModel(ref_model)

# Visualize Model Data
DimPlot(ref_data, group.by = "cell_type", label = TRUE,
        repel = TRUE) +
  ggtitle("scPred:: PBMC_1 Reference")

# Visualize predicted cell types
srat_scpred <- scPredict_edited(srat, ref_model)
DimPlot(srat_scpred, group.by = "scpred_prediction", label = TRUE,
        repel = TRUE) +
  ggtitle("Cell Types Predicted by scPred")




## Cell Classification with SingleR#############################################
#_______________________________________________________________________________
# Load dataset of immune cells bulk RNA-seq (platelets not included)
ref.se <- celldex::DatabaseImmuneCellExpressionData()
# Label celltypes in our srat dataset
pred.hesc <- SingleR::SingleR(test = srat@assays$RNA$counts, ref = ref.se,
                              assay.type.test=1,
                     labels = ref.se$label.fine)

# Add labels to metadata in srat object
srat@meta.data$singler_cell_types <- pred.hesc$pruned.labels
# UMAP Plot of Scitype annotated cells
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE,
        group.by = 'singler_cell_types') +
  ggtitle("SingleR with celldex::ImmuneDataset Ref")



# Exploratory Analysis (misc extra plots)#######################################
#_______________________________________________________________________________
# Visualize QC metrics as a violin plot
Idents(srat) <- "pbmc3k"
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,  idents = NULL, group.by = NULL,  split.by = NULL,   
        assay = "RNA")
# Compare QC features pairwise
plot1 <- FeatureScatter(srat, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(srat, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2


# Plot Variable features
plot1 <- VariableFeaturePlot(srat)
plot1
# Label Most variable features
plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(srat), 10), 
                     repel = TRUE)
plot2


# Visualize PCA Dim as scatter plot
VizDimLoadings(srat, dims = 1:2, reduction = "pca")

# Visualize PCA Dim as heatmap
DimHeatmap(srat, dims = 1:2, cells = 500, balanced = TRUE)

Idents(srat) <- srat$seurat_clusters
# Visualize gene expression of Top 2 Variable Genes across clusters
VlnPlot(srat, features = head(VariableFeatures(srat), 2))

# Visualize heatmap of genes across clusters
FeaturePlot(srat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                               "FCGR3A", "LYZ", "PPBP", "CD8A"))


# Heatmap of expression of top markers
topn <- srat_int.all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() 
DoHeatmap(srat, features = topn$gene) +
  NoLegend()



# DGE and Conserved Gene Expression ############################################
#_______________________________________________________________________________
# https://satijalab.org/seurat/archive/v3.1/immune_alignment.html
# We use the simulated integrated dataset we creat previously (randomly 
# assigning cells between (0) control group and (1) treatment group).
DefaultAssay(srat_int) <- "RNA"

# UMAP Plot of Scitype annotated cells
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE,
        group.by = 'cell_type') +
  ggtitle("SciType Annotated Cells")

# Visualize UMAP clusters
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE)


## Identify Conserved markers across conditions ################################
#_______________________________________________________________________________
conserved_marks <- FindConservedMarkers(srat_int, ident.1 = 1,   
                                          grouping.var = "group_id",
                                          verbose = FALSE)
head(conserved_marks)

### Visualize Top conserved markers for classical monocytes for all clusters 
#_______________________________________________________________________________
# Minimum cut-off set to 9th quantile
FeaturePlot(srat, features = rownames(head(conserved_marks)),
            min.cutoff = "q9")

### Visualize conserved marker expression with dot plot 
#_______________________________________________________________________________
DotPlot(srat, features = rev(rownames(conserved_marks[1:10,])), 
        cols = c("blue", "red"), dot.scale = 8,  split.by = "group_id") + 
  RotatedAxis()


## Differential Gene Expression: Option 1 (Naive) 
# Subset by each cell_type, find diff markers between conditions
#_______________________________________________________________________________
# Caution: With multiple samples, does not control for within sample variation
# Relabel cell identity label to cell_type (previously was cluster number)
cell_types <- levels(srat@meta.data$cell_type)
diff_markers = list()
for (n in seq_along(cell_types)) {
  # Isolate cells from first cell type/cluster
  sub_srat = subset(srat, idents = cell_types[n])
  # Reassign idents for finding markers
  Idents(sub_srat) = srat@meta.data$group_id
  diff_markers[[cell_types[n]]] <- 
    FindMarkers(sub_srat, ident.1 = "Ctrl", ident.2 = "Tx", slot = "scale.data")
  head(diff_markers[[cell_types[n]]], n = 10)
}

### Visualize diff marker expression with dot plot #############################
#_______________________________________________________________________________
Idents(srat) <- srat$cell_type
DotPlot(srat, features = rev(rownames(diff_markers$`Classical Monocytes`)[1:10]), 
        cols = c("blue", "red"), dot.scale = 8,  split.by = "group_id") + 
  RotatedAxis()

# Option 2: Visualize Differnetially Expressed Genes

# celltypes <- levels(Idents(srat))
# for (n in seq_along(celltypes)) {
#   sub_srat <- subset(srat, idents = celltypes[n])
#   Idents(sub_srat) <- "stim"
#   avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
#   avg.t.cells$gene <- rownames(avg.t.cells)
#   
# }


# Heatmap of gene expression between study groups across all cell types 
FeaturePlot(srat, features = c("LYZ", "ISG15"), 
            split.by = "group_id", max.cutoff = 3, 
            cols = c("grey", "red"))


# Visualize expression between study groups across all cell types
plots <- VlnPlot(srat, features = c("LYZ", "ISG15"), split.by = "group_id", 
                 group.by = "cell_type", pt.size = 0, combine = FALSE, 
                 split.plot = FALSE)
wrap_plots(plots = plots, ncol = 1)




## Differential Gene Expression: Option 3, Psuedo-bulk analysis ################
#_______________________________________________________________________________
# https://satijalab.org/seurat/articles/de_vignette

# Note: only works if tissue acquired from multiple replicates (not the case
#  with this dataset).So we simulate replicates.
srat$sample_id <- sample(x = 1:10, size = ncol(srat), replace = TRUE)

# Perform pseudo bulk, grouping gene expression by cell_type, group_id, 
#   and donor (can also group by sample/ donor if that exists in dataset)
pseudo_srat <- AggregateExpression(
  object = srat, assays = "RNA", return.seurat = T,
  group.by = c("group_id",  "cell_type", "sample_id"))
# For pseudo bulk testing we need to group by cell type and study group
pseudo_srat$celltype.tx <- paste(pseudo_srat$cell_type, 
                                 pseudo_srat$group_id, sep = "_")

# Set primrary identify/groups for cells for DGE test
Idents(pseudo_srat) <- "celltype.tx"
bulk.mono.de <- FindMarkers(object = pseudo_srat, 
                            ident.1 = "Classical Monocytes_Ctrl", 
                            ident.2 = "Classical Monocytes_Tx",
                            test.use = "DESeq2")
head(bulk.mono.de, n = 10)


### Visualize differentially expressed markers from pseudobulk analysis ########
#_______________________________________________________________________________
Idents(srat) <- srat$cell_type
DotPlot(srat, features = rev(rownames(bulk.mono.de)[1:10]), 
        cols = c("blue", "red"), dot.scale = 8,  split.by = "group_id") + 
  RotatedAxis()




# Monocle3 Pseudotime, tutorial dataset ########################################
#_______________________________________________________________________________

# Example small dataset (real data)
celegans_path <- here::here("_temp_data", "celegans_embryo", 
                            "count_matrix.Rdata")
if (!file.exists(celegans_path)) {
  expression_matrix <- 
    readRDS(url(paste0("https://depts.washington.edu:/trapnell-lab/software/",
                       "monocle3/celegans/data/packer_embryo_expression.rds")))
  cell_metadata <- 
    readRDS(url(paste0("https://depts.washington.edu:/trapnell-lab/software/",
                       "monocle3/celegans/data/packer_embryo_colData.rds")))
  gene_annotation <- 
    readRDS(url(paste0("https://depts.washington.edu:/trapnell-lab/software/",
                       "monocle3/celegans/data/packer_embryo_rowData.rds")))
  
  # Save expression data locally
  dir.create(here::here("_temp_data", "celegans_embryo"))
  save(expression_matrix, cell_metadata, gene_annotation,
       file = celegans_path)
}
load(file = celegans_path)


# Load 10x dataset into Monocle
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Set Random number seed for monocle
set.seed(1) 

# 1) Standard data normalization:
# A. Data Normalized by log and size factor for relative gene expression
# B. Performs linear dimension reduction
cds <- preprocess_cds(cds, num_dim = 100, method = "PCA")
# Plot variance explained with PCA
plot_pc_variance_explained(cds)
# 1B) Optional: scales datasets between runs/ timepoints (dataset integration)
cds <- 
  align_cds(cds, alignment_group = "batch", residual_model_formula_str =
              paste0("~ bg.300.loading + bg.400.loading + bg.500.1.loading",
                     " + bg.500.2.loading + bg.r17.loading + bg.b01.loading",
                     " + bg.b02.loading"))
# 2) Nonlinear dimension reduction (UMAP or tSNE) based on pre-processing
# umap.fast_sgd = FALSE: Makes the results reproducible
cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.fast_sgd = FALSE,
                        preprocess_method = 'PCA', cores = 1)


# 2A) Inspection: Visualize cells project onto reduced dimensional space
plot_cells(cds, label_groups_by_cluster = FALSE,  color_cells_by = "cell.type",
           show_trajectory_graph = FALSE)

# 2B) Inspection: Plot relative gene expression.
plot_cells(cds, genes=c("che-1", "hlh-17", "nhr-6", "dmd-6","ceh-36", "ham-1"), 
           label_cell_groups = FALSE, show_trajectory_graph = FALSE)

# 3) Assigns cells to clusters based on nonlinear dimension reduction projection
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 3A) Monocle needs cells groups by partition (metaclusters) as well as clusters
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE, 
           group_label_size = 4)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE,
           group_label_size = 4)

# 4) Fit principal graph (trajectory) for each partition
# Fits trajectories between projected cells
# Identifies branchpoints and endpoints with cell trajectories
cds <- learn_graph(cds)
# 4A) Visualize trajectories
plot_cells(cds, color_cells_by = "cell.type", label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, 
           group_label_size = 4, alpha = 0.5)

# 4B) Order cells by pseudotime and label branch points
plot_cells(cds, color_cells_by = "embryo.time.bin", label_cell_groups = FALSE,
           label_leaves = TRUE, label_branch_points=TRUE, graph_label_size = 3, 
           alpha = 0.5)
# 4C) Color cells by celltype and and label principle points (for ordering)
plot_cells(cds, color_cells_by = "cell.type", label_cell_groups = FALSE,
           label_leaves = TRUE, label_principal_points = TRUE, 
           graph_label_size = 3, alpha = 0.5)


# 5) Order cells in pseudotime from selected principle node(s)
cds <- order_cells(cds, root_pr_nodes = c("Y_261"))
# 5A) visualize pseudotime from root node
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, 
           graph_label_size = 1.5, alpha = 0.5)

# 6) Subset cells from a particular trajectory (graph segment)
cds_sub <- choose_graph_segments(cds, reduction_method = "UMAP",
                                 starting_pr_node = "Y_261", 
                                 ending_pr_nodes = c("Y_319"),
                                 clear_cds = TRUE) 

# Must repeat processing pipeline on this segment for next analysis
# 1) normalization and linear dimension reduction
cds_sub <- preprocess_cds(cds_sub, num_dim = 100, method = "PCA")
# 2) Nonlinear dimension reduction (UMAP) based on preprocessing
cds_sub <- reduce_dimension(cds_sub, reduction_method = "UMAP", 
                            preprocess_method = "PCA")
# 3) Cluster cells in reduced space
cds_sub <- cluster_cells(cds_sub)
# 4) Construct new graph
cds_sub <- learn_graph(cds_sub)
# 4A) Color cells by celltype and and label principle points
plot_cells(cds_sub, color_cells_by = "cell.type", 
           label_cell_groups = FALSE, label_leaves = TRUE, 
           label_principal_points = TRUE, alpha = 0.5,
           graph_label_size = 3)



# 5) Genes that are associated with the selected trajectory 
#  (cores = 1) for reproducibility
subset_pr_test_res <- graph_test(cds_sub, neighbor_graph="principal_graph", 
                                 cores = 1)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
# 6) Group genes into modules to visualize expression trends over pseudotime
gene_module_df <- find_gene_modules(cds_sub[pr_deg_ids,], resolution = 0.001)

# 7) Order modules by similarity (via hclust) to see which ones activate earlier
agg_mat <- aggregate_gene_expression(cds_sub, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)
                                [module_dendro$order])

# Visualize gene module activation over pseudotime
plot_cells(cds_sub, genes=gene_module_df, label_cell_groups = TRUE, 
           show_trajectory_graph = TRUE)





## Monocle with PBMC3K Dataset, Seurat Bridged #################################
#_______________________________________________________________________________

# Convert seurat object to monocle with wrapper (srat object from beginning)
cds <- SeuratWrappers::as.cell_data_set(srat)
# Bugfix from thread: https://github.com/satijalab/seurat-wrappers/issues/54
## Calculate size factors using built-in function in monocle3, add gene names
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- 
  rownames(srat[["RNA"]])

# Nonlinear dimension reduction
cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.fast_sgd = FALSE,
                        preprocess_method = 'PCA', cores = 1)

# Monocle needs partitions as well as clusters
# Using cluster_method = leiden raises error with Nonsymmetric adjacency matrix
# cluster_method = c("leiden", "louvain")
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 20, num_iter = 1,
                     cluster_method = "louvain", partition_qval = 0.05,
                     weight = FALSE, random_seed = 1,  verbose = FALSE)

# 3A) Monocle needs cells groups by partition (metaclusters) as well as clusters
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE, 
           group_label_size = 4)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE,
           group_label_size = 4)
plot_cells(cds, color_cells_by = "cell_type", show_trajectory_graph = FALSE,
           group_label_size = 4, alpha = 0.5)


# 4) Fit principal graph (trajectory) for each partition
# Fits trajectories between projected cells
# Identifies branchpoints and endpoints with cell trajectories
cds <- learn_graph(cds)
# 4A) Visualize trajectories
plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, 
           group_label_size = 4, alpha = 1, label_cell_groups = FALSE)

# # 4B) Order cells in pseudotime and label branch points
# plot_cells(cds, color_cells_by = "embryo.time.bin", label_cell_groups = FALSE,
#            label_leaves = TRUE, label_branch_points=TRUE, graph_label_size = 3, 
#            alpha = 0.5)
# 4C) Color cells in celltype and and label principle points
plot_cells(cds, color_cells_by = "cell_type", label_cell_groups = FALSE,
           label_leaves = TRUE, label_principal_points = TRUE, 
           graph_label_size = 3, alpha = 0.5)

# NOTE: The node ID will chang between systems, you need to look at the graph 
# above and change the node id "Y_xxx" and "Y_xxx" in the function calls below.

# 5) Order cells in pseudotime from selected principle node(s)
cds <- order_cells(cds, root_pr_nodes = c("Y_311"))
# 5A) visualize pseudotime from root node
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE, 
           graph_label_size = 1.5, alpha = 0.5)

# # 6) Subset cells from a particular trajectory (graph segment)
# cds_sub <- choose_graph_segments(cds, reduction_method = "UMAP",
#                                  starting_pr_node = "Y_311", 
#                                  ending_pr_nodes = c("Y_324"),
#                                  clear_cds = TRUE) 
# 
# # Must repeat processing pipeline on this segment for next analysis
# # 1) normalization and linear dimension reduction
# cds_sub <- preprocess_cds(cds_sub, num_dim = 100, method = "PCA")
# # 2) Nonlinear dimension reduction (UMAP) based on preprocessing
# cds_sub <- reduce_dimension(cds_sub, reduction_method = "UMAP", 
#                             preprocess_method = "PCA")
# # 3) Cluster cells in reduced space
# cds_sub <- cluster_cells(cds_sub)
# # 4) Construct new graph
# cds_sub <- learn_graph(cds_sub)
# # 4A) Color cells in celltype and and label principle points
# plot_cells(cds_sub, color_cells_by = "cell_type", 
#            label_cell_groups = FALSE, label_leaves = TRUE, 
#            label_principal_points = TRUE, alpha = 0.5,
#            graph_label_size = 3)
# 
# 
# 
# # 5) Genes that are associated with the selected trajectory 
# #  (cores = 1) for reproducibility
# subset_pr_test_res <- graph_test(cds_sub, neighbor_graph="principal_graph", 
#                                  cores = 1)
# pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
# # 6) Group genes into modules to visualize expression trends over pseudotime
# gene_module_df <- find_gene_modules(cds_sub[pr_deg_ids,], resolution = 0.001)
# 
# # 7) Order modules by similarity (via hclust) to see which ones activate earlier
# agg_mat <- aggregate_gene_expression(cds_sub, gene_module_df)
# module_dendro <- hclust(dist(agg_mat))
# gene_module_df$module <- factor(gene_module_df$module, 
#                                 levels = row.names(agg_mat)
#                                 [module_dendro$order])
# 
# # Visualize gene module activation over pseudotime
# plot_cells(cds_sub, genes=gene_module_df, label_cell_groups = TRUE, 
#            show_trajectory_graph = TRUE)
# 
