


library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(HGNChelper)
library(openxlsx)
library(ggplot2)

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
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = here("_temp_data", "filtered_gene_bc_matrices/hg19"))

# Initialize the Seurat object with the raw count matrix (non-normalized data).
# Threshold genes that are found within at least 3 cells
# Threshold cells that have at least 200 genes
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Add column in metadata slot for percent of mitochondrial genes (QC metric)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# nFeature_RNA > 200: removes empty droplets or poor quality cells with little DNA
# nFeature_RNA <25000: remove doublets (droplets with 2+ cells)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
# 1. Normalizes gene expression by the total expression in each cell
# 2. Multiplies this by a scale factor (10,000 by default)
# 3. Log-transforms the result.
# Stored in: pbmc[["RNA"]]$data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature Selection
# Identify highly variables genes, use these as anchors for downstream dimension
# reduction
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data (across cells, only selected variable features)
# Essentially converts gene expression to z-score (normalize by mean and std across cells)
# Stored in: pbmc[["RNA"]]$scale.data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# Scale data can also be used to remove unwanted cell cycle variation
# However, this is a more advance method, and it is recommended to use the new
# Seurat workflow: SCTransform(). 
# Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9
# Vignette: https://satijalab.org/seurat/articles/sctransform_vignette
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")


# Linear dimension reduction (PCA)
#-------------------------------------------------------------------------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Plot commands: VizDimReduction(), DimPlot(), and DimHeatmap()
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca").     ** most useful**
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# DimPlot(pbmc, reduction = "pca") + NoLegend()

# Choose dimensionality of the dataset
# Maximize the signal (biological variability) to the noise (other sources of variation)
ElbowPlot(pbmc)

# Clustering
#-------------------------------------------------------------------------------
# Construct a kNN graph based on euclidean distance in a subset of PCA space 
#  (up to dimensionality chosen).
# Refine edge weights between pairs of cells based on their shared overlap and 
# local neighboors (Jaccard similarity)
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# To cluster the cells, we next apply modularity optimization (Louvain algorithm, SLM )
#  to iteratively group cells together, with the goal of optimizing the standard 
#  modularity function. 
# Cluster granularity is set with resolution, 0.4-1.2 typically returns good results
#   for single-cell datasets of around 3K cells. 
#   resolution often increases for larger datasets.
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Perform UMAP clustering
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize UMAP clusters
DimPlot(pbmc, reduction = "umap")

# Cluster Marker Identification
#-------------------------------------------------------------------------------
# Note: Install presto package for faster results
# Findconservedmarkers()
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)




# Scitype Cell Type Identification
#-------------------------------------------------------------------------------
#https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md
# load gene set and cell type annotation functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" 
# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,
# Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = pbmc[["RNA"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either 
# 1. pbmc[["RNA"]]@scale.data (default), 
# 2. pbmc[["SCT"]]@scale.data, if sctransform is used for normalization,
# 3. pbmc[["integrated"]]@scale.data, for joint analysis of multiple sc datasets.

# Merge by cluster
cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# Set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# Add another column in seurat metadata for celltype
pbmc@meta.data$cell_type_scitype <- select(pbmc@meta.data, "seurat_clusters") %>%
  left_join(y = select(sctype_scores, "cluster", "type"), 
            by = join_by(seurat_clusters == cluster)) %>% pull("type")

# UMAP Plot of Scitype annotated cells
DimPlot(pbmc, reduction = "umap", label = FALSE, repel = TRUE, 
        group.by = 'cell_type_scitype') + 
  ggtitle("SciType Annotated Cells")      


# Cell Trajectory Analysis


