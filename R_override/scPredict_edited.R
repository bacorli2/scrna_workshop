.make_names <- function(x){
  x <- gsub("\\+", "_plus", x)
  x <- gsub("\\-", "_minus", x)
  x <- make.names(x)
}

scPredict_edited <- function (new, reference, threshold = 0.55, max.iter.harmony = 20, 
                              recompute_alignment = TRUE, seed = 66) 
{
  if (!(is(reference, "Seurat") | is(reference, "scPred"))) 
    stop("'object' must be of class 'scPred' or 'Seurat'")
  if (is(reference, "Seurat")) {
    spmodel <- reference@misc$scPred
  }
  else {
    spmodel <- reference
  }
  if (is.null(spmodel)) 
    stop("No feature space has been determined!")
  if (!length(spmodel@train)) 
    stop("No models have been trained!")
  if (!is(new, "Seurat")) 
    stop("New data must be a Seurat object")
  new <- project_query_edited(new, reference = spmodel, max.iter.harmony = max.iter.harmony, 
                              recompute_alignment = recompute_alignment, seed = seed)
  new_embeddings_aligned <- Embeddings(new[["scpred"]])
  colnames(new_embeddings_aligned) <- colnames(spmodel@cell_embeddings)
  cellTypeModelNames <- names(spmodel@features)
  .predictCellClass <- function(cellType, spmodel, testEmbeddings) {
    features <- as.character(spmodel@features[[cellType]]$feature)
    model <- spmodel@train[[cellType]]
    prediction <- predict(model, newdata = scPred:::subsetMatrix(testEmbeddings, 
                                                                 features), type = "prob")
    rownames(prediction) <- rownames(testEmbeddings)
    prediction[, 1, drop = FALSE]
  }
  cat(crayon::green(cli::symbol$record, " Classifying cells...\n"))
  res <- sapply(cellTypeModelNames, .predictCellClass, spmodel, 
                new_embeddings_aligned)
  res <- as.data.frame(res)
  colnames(res) <- cellTypeModelNames
  rownames(res) <- colnames(new)
  classes <- cellTypeModelNames
  if (length(cellTypeModelNames) == 1) {
    metadata <- get_metadata(spmodel)
    cellClasses <- levels(metadata$pvar)
    res_comp <- 1 - res[, 1]
    negClass <- cellClasses[cellClasses != names(res)]
    res[[negClass]] <- res_comp
  }
  max_props <- as.data.frame(t(apply(res, 1, function(x) c(index = which.max(x), 
                                                           max = x[which.max(x)]))))
  names(max_props) <- c("index", "max")
  max_props$generic_class <- names(res)[max_props$index]
  res <- cbind(res, max_props)
  pred <- ifelse(res$max > threshold, res$generic_class, "unassigned")
  names(pred) <- colnames(new)
  res$prediction <- pred
  res$index <- NULL
  res$no_rejection <- res$generic_class
  res$generic_class <- NULL
  names(res) <- .make_names(paste0("scpred_", names(res)))
  new <- AddMetaData(new, res)
  cat(crayon::green("DONE!\n"))
  new
}

project_query_edited <- function (new, reference, max.iter.harmony = 20, recompute_alignment = TRUE, 
                                  seed = 66, ...) 
{
  if (!(is(reference, "Seurat") | is(reference, "scPred"))) 
    stop("'object' must be of class 'scPred' or 'Seurat'")
  if (is(reference, "Seurat")) {
    spmodel <- reference@misc$scPred
  }
  else {
    spmodel <- reference
  }
  if (is.null(spmodel)) 
    stop("No feature space has been determined!")
  if (!is(new, "Seurat")) 
    stop("New data must be a Seurat object")
  if ("scpred" %in% names(new@reductions)) {
    if (recompute_alignment) {
      alignment <- TRUE
      cat(crayon::yellow(cli::symbol$figure_dash, "Data has already been aligned to a reference.\n"), 
          sep = "")
      cat(crayon::yellow(cli::symbol$sup_plus, "Skip data alignment using `recompute.alignment = FALSE`.\n"), 
          sep = "")
    }
    else {
      alignment <- FALSE
    }
  }
  else {
    alignment <- TRUE
  }
  if (alignment) {
    cat(crayon::green(cli::symbol$record, " Matching reference with new dataset...\n"))
    ref_loadings <- spmodel@feature_loadings
    ref_embeddings <- spmodel@cell_embeddings
    new_features <- rownames(new)
    reference_features <- rownames(ref_loadings)
    shared_features <- intersect(reference_features, new_features)
    cat(crayon::cyan("\t", cli::symbol$line, paste(length(reference_features), 
                                                   "features present in reference loadings\n")))
    cat(crayon::cyan("\t", cli::symbol$line, paste(length(shared_features), 
                                                   "features shared between reference and new dataset\n")))
    cat(crayon::cyan("\t", cli::symbol$line, paste0(round(length(shared_features)/length(reference_features) * 
                                                            100, 2), "% of features in the reference are present in new dataset\n")))
    ref_loadings <- ref_loadings[shared_features, ]
    new_data <- GetAssayData(new, layer = "data")[shared_features, 
    ]
    means <- spmodel@scaling$means
    stdevs <- spmodel@scaling$stdevs
    new_data <- Matrix::t(new_data)
    names(means) <- names(stdevs) <- rownames(spmodel@scaling)
    means <- means[shared_features]
    stdevs <- stdevs[shared_features]
    i <- stdevs == 0
    if (any(i)) {
      warning(paste0(sum(i), " features have zero variance but are present in the feature loadings. \nDid you subset or integrated this data before?"))
      cat(crayon::yellow("Removing zero-variance genes from projection\n"))
      new_data <- new_data[, !i]
      ref_loadings <- ref_loadings[!i, ]
      means <- means[!i]
      stdevs <- stdevs[!i]
    }
    scaled_data <- scale(new_data, means, stdevs)
    new_embeddings <- scaled_data %*% ref_loadings
    dataset <- factor(c(rep("reference", nrow(ref_embeddings)), 
                        rep("new", nrow(new_embeddings))), levels = c("reference", 
                                                                      "new"))
    rownames(ref_embeddings) <- paste0("ref_", rownames(ref_embeddings))
    rownames(new_embeddings) <- paste0("new_", rownames(new_embeddings))
    eigenspace <- as.data.frame(rbind(ref_embeddings, new_embeddings))
    meta_data <- data.frame(rownames(eigenspace), dataset = dataset)
    cat(crayon::green(cli::symbol$record, " Aligning new data to reference...\n"))
    set.seed(seed)
    harmony_embeddings <- harmony::HarmonyMatrix(eigenspace, meta_data, 
                                                 "dataset", do_pca = FALSE, reference_values = "reference", 
                                                 max.iter.harmony = max.iter.harmony, ...)
    new_embeddings_aligned <- harmony_embeddings[dataset == 
                                                   "new", , drop = FALSE]
  }
  else {
    new_embeddings_aligned <- Embeddings(new, reduction = "scpred")
    colnames(new_embeddings_aligned) <- gsub("scpred_", spmodel@reduction_key, 
                                             colnames(new_embeddings_aligned))
  }
  rownames(new_embeddings_aligned) <- gsub("^new_", "", rownames(new_embeddings_aligned))
  new@reductions[["scpred"]] <- CreateDimReducObject(embeddings = new_embeddings_aligned, 
                                                     key = "scpred_", assay = DefaultAssay(object = new))
  if (recompute_alignment) {
    rownames(new_embeddings) <- gsub("^new_", "", rownames(new_embeddings))
    new@reductions[["scpred_projection"]] <- CreateDimReducObject(embeddings = new_embeddings, 
                                                                  key = "Projection_", assay = DefaultAssay(object = new))
  }
  new
}
