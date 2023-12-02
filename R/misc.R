#' Calculates QPSP profiles of samples
#' Rownames of df1 has to be the same annotation type as list of network vectors
#' @param df1 dataframe of discrete GFS-transformed expression data with
#' features in rows and samples in columns
#' @param list_complex list of vectors that contain components of networks
#' @return dataframe of QPSP profiles
calcQPSP <- function(df1, list_complex) {
  # All components in complex have to be measured in df1
  df1_proteins <- rownames(df1)
  # Logvec of whether entire complex is present
  complex_logvec <- sapply(list_complex,
                           function(vec) all(vec %in% df1_proteins))
  print(sum(complex_logvec))
  # Error if any NA in complex_logvec
  if (anyNA(complex_logvec)) stop("NA present in logvec")
  sublist_complex <- list_complex[complex_logvec]
  #' Helper function that calculates QPSP profile for single sample
  #' @param col_vec Column vector representing single sample
  calcSingleQPSP <- function(col_vec) {
    sapply(sublist_complex, function(vec) mean(col_vec[vec]))
  }
  return(data.frame(apply(df1, 2, calcSingleQPSP)))
}


# Evaluates DE analysis method
# Arguments: 2 logical vectors - predicted and truth labels
# Returns: Accuracy, Sensitivity, Precision, Specificity, NPR
evaluation_report <- function(predict_vec, label_vec) {
  accuracy <- sum(predict_vec == label_vec)/length(predict_vec)
  TP <- sum(predict_vec & label_vec)
  TN <- sum(!predict_vec & !label_vec)
  sensitivity_1 <- TP/sum(label_vec)
  precision_1 <- TP/sum(predict_vec)
  specificity_1 <- TN/sum(!label_vec)
  negative_predictive_rate_1 <- TN/sum(!predict_vec)
  Sensitivity <- c(sensitivity_1, specificity_1)
  Precision <- c(precision_1, negative_predictive_rate_1)
  df <- cbind(Sensitivity, Precision)
  rownames(df) <- c("DE", "Not DE")
  print(df)
  print(paste("Accuracy:", accuracy))
  metrics <- c(accuracy, sensitivity_1, precision_1,
               specificity_1, negative_predictive_rate_1)
  names(metrics) <- c("Accuracy", "Sensitivity", "Precision",
                      "Specificity", "NPR")
  return(metrics)
}


# Side effects: Writes KEGG pathway dataframes (Entrez)
kegg_df <- function(kegg_fpath) {
  pathway_id <- substring(kegg_fpath, 29, 36)
  wpath <- sprintf("../info/KEGG/hsa_df/%s.tsv", pathway_id)
  # Parses file into dataframe
  df_hsa <- parseKGML2DataFrame(kegg_fpath)
  # Convert KEGG ID to Entrez ID
  df_hsa$from <- substring(df_hsa$from, 5)
  df_hsa$to <- substring(df_hsa$to, 5)
  write.table(df_hsa, wpath,
              quote = F, sep = "\t", row.names = F)
}

#' Evaluate batch effects using a given set of metrics
#'
#' @param obj Either a Seurat object or dataset obj 
eval_batch <- function(
  obj, batchname, classname, slot = "data",
  metrics = c("rvp", "cms", "kbet", "lisi", "gpca", "pvca"),
  nperm = 0, k.cms = 50, k0 = NULL, perplexity = 30,
  ret.scores = FALSE
) {
  rvp_obj <- sce_coldata <- kbet <- lisi <- gpca <- pvca <- NULL
  scores <- rep(NA, 6)
  names(scores) <- c("rvp", "cms", "kbet", "lisi", "gpca", "pvca")

  if (class(obj) == "Seurat") {
    cat("Detected: Seurat object", fill = TRUE)
    X <- Seurat::GetAssayData(obj, slot = slot)
    metadata <- obj@meta.data
    n <- dim(obj)[2]
  } else if (class(obj) == "SingleCellExperiment") {
    cat("Detected: SingleCellExperiment object", fill = TRUE)
    X <- SingleCellExperiment::logcounts(obj)
    metadata <- obj@colData
    n <- ncol(obj)
  } else {
    cat("Detected: Dataset object", fill = TRUE)
    X_raw <- obj$X
    X <- log2_transform(scale_trimmed(2 ^ X_raw)) # log-normalise
    metadata <- obj$metadata
    n <- ncol(X)
  }

  for (metric in metrics) {
    if (metric == "rvp") {
      cat("Calculating RVP...", fill = TRUE)
      rvp_obj <- rvp(
        as.array(Matrix::t(X)),
        metadata[[batchname]],
        metadata[[classname]],
        ret.percent = FALSE
      )
      scores[metric] <- rvp_obj$percent.batch
    } else if (metric == "cms") {
      cat("Calculating CMS...", fill = TRUE)
      library(CellMixS)
      sce <- SingleCellExperiment::SingleCellExperiment(
        list(logcounts = X), colData = metadata
      )
      sce <- cms(sce, k = k.cms, group = batchname)
      sce_coldata <- sce@colData
      scores[metric] <- mean(sce@colData$cms)
    } else if (metric == "kbet") {
      cat("Calculating kBET...", fill = TRUE)
      library(kBET)
      kbet <- kBET(
        Matrix::t(X),
        metadata[[batchname]],
        k0 = k0,
        testSize = n,
        n_repeat = 1,
        verbose = TRUE
      )
      cat("\n")
      scores[metric] <- kbet$summary$kBET.observed[1]
    } else if (metric == "lisi") {
      cat("Calculating LISI...", fill = TRUE)
      library(lisi)
      lisi <- compute_lisi(
        Matrix::t(X),
        metadata,
        c(batchname, classname),
        perplexity = perplexity
      )
      scores[metric] <- mean(lisi[[batchname]])
    } else if (metric == "gpca") {
      cat("Calculating gPCA...", fill = TRUE)
      # modified to fix error when nperm = 0
      gpca <- gPCA(Matrix::t(X), metadata[[batchname]], nperm = 0)
    } else if (metric == "pvca") {
      cat("Calculating PVCA...", fill = TRUE)
      library(pvca)
      meta_metadata <- data.frame(labelDescription = colnames(metadata))
      pheno_data <- new(
        "AnnotatedDataFrame", data = metadata, varMetadata = meta_metadata
      )
      eset <- Biobase::ExpressionSet(assayData = X, phenoData = pheno_data)
      pvca <- pvcaBatchAssess(eset, c('batch', 'class'), 0.6)
    }
  }

  if (ret.scores) {
    return(scores)
  } else {
    return(list(
      rvp = rvp_obj,
      cms = sce_coldata,
      kbet = kbet,
      lisi = lisi,
      gpca = gpca,
      pvca = pvca
    ))
  }
}


eval_batch_effects <- function(df, batch_info, class_info) {
  # Tranpose df to n x p
  pca_obj <- prcomp(t(df), center = T, scale. = F)
  # Eigenvalues
  pca_df <- data.frame(pca_obj$x)
  eig_value <- (pca_obj$sdev)^2
  # Percentage variance of each PC
  var_pct <- eig_value/sum(eig_value)
  
  # Libraries: cluster, gpca
  # Argument: PC coordinates for single PC
  calc_pc_var <- function(vec) {
    # Argument: ANOVA attributes; Calculates percentage of variance
    # SS_between/(SS_between + SS_within)
    calc_var_percentage <- function(vec) unname(vec[3]/(vec[3] + vec[4]))
    pc_metadata <- data.frame(pc = vec,
                              batch_info = as.factor(batch_info),
                              class_info = as.factor(class_info))
    batch_anova_attr <- unlist(summary(aov(pc~batch, data = pc_metadata)))
    class_anova_attr <- unlist(summary(aov(pc~class, data = pc_metadata)))
    return(c(calc_var_percentage(batch_anova_attr),
             calc_var_percentage(class_anova_attr)))
  }
  var_composition <- sapply(pca_df, calc.pc_var)
  # Dataframe showing prop. var and pi_BE and pi_BV
  pca_var_df <- data.frame(var_pct,
                           batch_pct = var_composition[1,],
                           class_pct = var_composition[2,])
  total_batch_pct <- sum(pca_var_df$var_pct * pca_var_df$batch_pct)
  total_class_pct <- sum(pca_var_df$var_pct * pca_var_df$class_pct)
  print(head(pca_var_df))
  
  # Calculating average silhouette scores
  dist_pca <- dist(pca_df, "euclidean")
  batch_df <- cluster::silhouette(batch_info, dist_pca)
  class_df <- cluster::silhouette(class_info, dist_pca)
  
  # gPCA
  # gpca_obj <- gPCA.batchdetect(t(df), batch_info)
  
  # # Collating results
  # metrics <- c(total_batch_pct*100, total_class_pct*100,
  #              mean(batch_df[,3]), mean(class_df[,3]),
  #              gpca_obj$delta, gpca_obj$p.val)
  # Collating results
  metrics <- c(total_batch_pct*100, total_class_pct*100,
               mean(batch_df[,3]), mean(class_df[,3]))
  
  # names(metrics) <- c("pi_BE", "pi_BV",
  #                     "silhouette_BE", "silhouette_BV",
  #                     "gPCA_delta", "gPCA_pvalue")[]
  print(xtable(t(matrix(metrics)), digits = 3))
  cat(metrics, sep = "\t")
  return(metrics)
}


calc_var_preservation <- function(df_bef, df_aft) {
  return(sum(apply(df_aft, 1, var))/sum(apply(df_bef, 1, var)))
}


#' Calculates sum of squares in a vector
#'
#' @param x numeric vector
#' @param x data.frame or matrix with dim (n_features, n_samples) 
sum_squares <- function(x) {
  if (is.vector(x)) {
    return(sum((x - mean(x)) ^ 2))
  } else if (is.matrix(x) | is.data.frame(x)) {
    feat_means <- rowMeans(x)
    squares <- (x - feat_means) ^ 2
    return(sum(squares))
  }
}


#' Calculates proportion of variance in dataframe due to batch effects and biological variable of interest
#' 
#' @param df p x n gene expression dataframe (p: no. of probes, n: no. of samples)
#' @param batch_info vector containing batch labels of samples (ordered the same way as in the dataframe)
#' @param class_info vector containing class labels of samples (ordered the same way as in the dataframe)
#' @return vector containing total proportion of variance in dataframe due to batch effects and biological variable of interest, respectively
calc_var_prop <- function(
  df, metadata, batch_name = "batch_info", class_name = "class_info"
) {
  # Tranpose df to n x p
  pca_obj <- prcomp(t(df))
  # Eigenvalues
  pca_df <- data.frame(pca_obj$x)
  eig_value <- (pca_obj$sdev)^2
  # Percentage variance of each PC
  var_pct <- eig_value/sum(eig_value)
  
  # Libraries: cluster, gpca
  # Argument: PC coordinates for single PC
  calc_pc_var <- function(vec) {
    # Argument: ANOVA attributes; Calculates percentage of variance
    # SS_between/SS_between + SS_within
    calc_var_percentage <- function(vec) unname(vec[3]/(vec[3] + vec[4]))
    pc_metadata <- data.frame(
      pc = vec,
      batch = metadata[rownames(pca_df), batch_name],
      class = metadata[rownames(pca_df), class_name]
    )
    batch_anova_attr <- unlist(summary(aov(pc~batch, data = pc_metadata)))
    class_anova_attr <- unlist(summary(aov(pc~class, data = pc_metadata)))
    return(c(calc_var_percentage(batch_anova_attr),
             calc_var_percentage(class_anova_attr)))
  }
  var_composition <- sapply(pca_df, calc_pc_var)
  # Dataframe showing prop. var and pi_BE and pi_BV
  pca_var_df <- data.frame(var_pct,
                           batch_pct = var_composition[1,],
                           class_pct = var_composition[2,])
  total_batch_pct <- sum(pca_var_df$var_pct * pca_var_df$batch_pct)
  total_class_pct <- sum(pca_var_df$var_pct * pca_var_df$class_pct)
  metrics <- c(total_batch_pct, total_class_pct)
  names(metrics) <- c("var_batch", "var_class")
  return(metrics)
}


#' Sums up distances between batches for all classes in PCA transformed data
#'
#' @param X dataframe of dim (n_features, n_samples)
calc_batch_dist <- function(
  X, metadata, pca = TRUE, batch_name = "batch_info", class_name = "class_info"
) {
  if (pca) {
    pca_obj <- prcomp(t(X))
    X_pca <- data.frame(pca_obj$x)
  } else {
    X_pca <- t(X)
  }
  dists <- numeric(ncol(X_pca))
  names(dists) <- colnames(X_pca)
  for (i in seq(ncol(X_pca))) {
    pc_i <- data.frame(
      pc = X_pca[, i],
      batch = metadata[rownames(X_pca), batch_name],
      class = metadata[rownames(X_pca), class_name]
    )
    classes <- unique(pc_i$class)
    dist_i <- 0
    for (c in classes) {
      pc_class <- pc_i[pc_i$class == c, ]
      batches <- unique(pc_class$batch)
      pc_batches <- lapply(batches, function(b) pc_class[pc_class$batch == b, ])
      medians <- sapply(pc_batches, function(X) median(X[, 1]))
      dist_mat <- dist(medians)^2
      dist_i <- dist_i + sum(dist_mat)
    }
    dists[i] <- dist_i 
  }

  dists
}
