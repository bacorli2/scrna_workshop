


library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(HGNChelper)
library(openxlsx)
library(presto)
library(ggplot2)
library(scAnnotatR)
library(SingleR)
library(celldex)

# This code was appropriate from the following sources
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

# Set wd to base of workshop repository
here::i_am("README.md")

# Example small dataset (real data)
# Used from this tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial
# 2,700 single cells that were sequenced on the Illumina NextSeq 500
# 13,714 genes
# Download dataset into temp_data, unzip
if (!file.exists(here::here("_temp_data", "pbmc3k_filtered_gene_bc_matrices.tar.gz"))) {
  dir.create(here::here("_temp_data"))
  download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                destfile = here::here("_temp_data", "pbmc3k_filtered_gene_bc_matrices.tar.gz"))
  untar(here::here("_temp_data", "pbmc3k_filtered_gene_bc_matrices.tar.gz"), 
        exdir = here::here("_temp_data"))
}
# Load the srat dataset
srat.data <- Read10X(data.dir = here("_temp_data", "filtered_gene_bc_matrices/hg19"))

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

# Filter poor quality cells
# nFeature_RNA > 200: removes empty droplets or cells with little RNA
# nFeature_RNA < 25000: remove doublets (droplets with 2+ cells)
# percent.mt < 5: removes cells with over 5% mitochondrial DNA 
# (poor viability)
srat <- subset(srat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                 percent.mt < 5)

# Normalize data
# 1. Normalizes gene expression by the total expression in each cell
# 2. Multiplies this by a scale factor (10,000 by default)
# 3. Log-transforms the result.
# Stored in: srat[["RNA"]]$data
srat <- NormalizeData(srat, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Feature Selection
# Identify highly variables genes, (to be used for dimension reduction)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)

# Scale the data (across cells, only selected variable features)
# Essentially converts gene expression to z-score (normalize by mean and std 
# across cells)
# Stored in: srat[["RNA"]]$scale.data
srat <- ScaleData(srat, features = rownames(srat))

# Scale data can also be used to remove unwanted cell cycle variation
# However, this is a more advance method, and it is recommended to use the new
# Seurat workflow: SCTransform(). 
# Paper: https://genomebiology.biomedcentral.com/articles/10.1186/
#        s13059-021-02584-9
# Vignette: https://satijalab.org/seurat/articles/sctransform_vignette
# srat <- ScaleData(srat, vars.to.regress = "percent.mt")


# Linear dimension reduction (PCA)
#-------------------------------------------------------------------------------
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
# Plot commands: VizDimReduction(), DimPlot(), and DimHeatmap()
VizDimLoadings(srat, dims = 1:2, reduction = "pca")     
DimHeatmap(srat, dims = 1, cells = 500, balanced = TRUE)
DimPlot(srat, reduction = "pca") + NoLegend()

# Choose dimensionality of the dataset
# Maximize the signal (biological variability) to the noise (other sources of 
# variation)
ElbowPlot(srat)


# Clustering
#-------------------------------------------------------------------------------
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



# Cluster Marker Identification
#-------------------------------------------------------------------------------

# Find differentially expressed genes in each cluster vs. all other clusters
# Test used is non-parametric Wilcoxon rank sum test
# Note: Install presto package for much faster results
srat.all.markers <- FindAllMarkers(srat, only.pos = TRUE)

# Findconservedmarkers()


# Cell Type Annotation: Scitype 
#-------------------------------------------------------------------------------
#https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md

# Load gene set and cell type annotation functions
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
srat@meta.data$cell_type_scitype <- select(srat@meta.data, "seurat_clusters") %>%
  left_join(y = select(sctype_scores, "cluster", "type"), 
            by = join_by(seurat_clusters == cluster)) %>% pull("type")

# UMAP Plot of Scitype annotated cells
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE, 
        group.by = 'cell_type_scitype') + 
  ggtitle("SciType Annotated Cells")      


# Cell Classification Using scAnnotateR
#-------------------------------------------------------------------------------

# DEFAULT MODEL: Load classification Models
default_models <- scAnnotatR::load_models("default")

# Perform classification
srat_scannot <- classify_cells(classify_obj = srat, 
                             assay = 'RNA', slot = 'counts',
                             cell_types = "all", 
                             path_to_models = 'default')
# Plot best match for each cell_type
DimPlot(srat_scannot, group.by = "most_probable_cell_type")


# Classification Based Cell Type: scPred
--------------------------------------------------------------------------------
 
source(here::here("R_override", "scPredict_edited.R"))
devtools::install_github("immunogenomics/harmony")
devtools::install_github("powellgenomicslab/scPred")

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




# Cell Classification with SingleR
#------------------------------------------------------------------------------
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




# Pseudotime with monocle3
#-------------------------------------------------------------------------------
# Code appropriated from:
# https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-
# in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
library(monocle3)

remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install("monocle", force = TRUE)

# BiocManager::install("GenomeInfoDb", force = TRUE)
# devtools::install_github("cysouw/qlcMatrix")

library(SeuratWrappers)


# Convert seurat object to monocle
cds <- SeuratWrappers::as.cell_data_set(srat)
# Monocle needs partitions as well as clusters
cds <- monocle3::cluster_cells(cds, resolution=1e-3)

# VIsualize clusters and larger partitions
p1 <- monocle3::plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- monocle3::plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
patchwork::wrap_plots(p1, p2)

# Subsetting partitions
integrated.sub <- base::subset(Seurat::as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- SeuratWrappers::as.cell_data_set(integrated.sub)

# Trajectory analysis
cds <-  monocle3::learn_graph(cds, use_partition = TRUE, verbose = FALSE)


# Visualize Trajectory
monocle3::plot_cells(cds, color_cells_by = "cluster",  label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)


# Color Cells by pseudo time
cds <-  monocle3::order_cells( cds, root_cells = colnames(  
  cds[, monocle3::clusters( cds) == 1]) )
monocle3::plot_cells(cds, color_cells_by = "pseudotime", group_cells_by = "cluster",
           label_cell_groups = FALSE, label_groups_by_cluster = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE,
           label_roots = FALSE, trajectory_graph_color = "grey60")

integrated.sub <- Seurat::as.Seurat(cds, assay = NULL)
# Monocle assigns cells outside of timeline and InF (erros with seurat)
# Reassign to NA to make seurat compatible
integrated.sub@meta.data$monocle3_pseudotime[is.infinite(integrated.sub@meta.data$monocle3_pseudotime)] <- NA
Seurat::FeaturePlot(integrated.sub, "monocle3_pseudotime")

# Detect genes that vary over a trajectory
# Can only use multicore on mac or linux
cds_graph_test_results <- 
  monocle3::graph_test(cds, neighbor_graph = "principal_graph",
                       cores = 1)
# If rbind error:
# trace(‘calculateLW’, edit = T, where = asNamespace(“monocle3”))
# find Matrix::rBind and replace with rbind then save.



# SummarizedExperiment::rowData(cds)$gene_short_name <- SummarizedExperiment::rowData(cds)$gene_name

# VIsualize top genes that varied over trajectorty
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)


deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.2))


# VIsualize most significant genes alog trajectorty
# Erroring
monocle3::plot_cells(cds,
           genes=head(deg_ids, n=1),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)


# "IFNG" %in% rownames(SummarizedExperiment::rowData(cds))   # TRUE
# "GZMB" %in% rownames(SummarizedExperiment::rowData(cds))    # TRUE
# 
# "IFNG" %in% SummarizedExperiment::rowData(cds)$gene_name    # TRUE
# "GZMB" %in% SummarizedExperiment::rowData(cds)$gene_name    # TRUE


library(monocle3)
gene_modules <- monocle3::find_gene_modules(cds[deg_ids,],
                                  resolution=c(10^seq(-6,-1)))
table(gene_modules$module)





