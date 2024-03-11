

scType_SeuratObj <- function (srat, tissue, cell_type_col = "Cell_type") {
  
  ## Cell Type Annotation: ScType ###############################################
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
  Idents(srat) <- srat@meta.data[[cell_type_col]]
  
  # UMAP Plot of Scitype annotated cells
  DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE,
          group.by = cell_type_col) +
    ggtitle("SciType Annotated Cells")
  
  return(srat)
  
}