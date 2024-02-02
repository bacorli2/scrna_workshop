


library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(HGNChelper)
library(openxlsx)
library(presto)
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


function(cl){
  es.max.cl = sort(rowSums(
    es.max[ ,rownames(srat@meta.data[srat@meta.data$seurat_clusters==cl, ])]),
    decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, 
                  ncells = sum(srat@meta.data$seurat_clusters==cl)), 10)
}

filter()

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
# Sctype scores scale by n, so threshold is score/4
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < 
                     sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# Add column in seurat metadata for celltype annotation
srat@meta.data$ctype_scitype <- select(srat@meta.data, "seurat_clusters") %>%
  left_join(y = select(sctype_scores, "cluster", "type"), 
            by = join_by(seurat_clusters == cluster)) %>% pull("type")

# UMAP Plot of Scitype annotated cells
DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE, 
        group.by = 'ctype_scitype') + 
  ggtitle("SciType Annotated Cells")      


# Cell Trajectory Analysis


