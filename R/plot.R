#' Plots PCA plot using ggplot2
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom Matrix t
#' @param X Data frame with dim (n_features, n_samples) or data frame of PC
#'   values with dim (n_samples, n_features)
#' @param metadata Data frame of metadata with dim (n_samples, n_features)
#' @param var_pc Numeric vector of variance explained by each PC 
#' @param newdata Dataframe of data to be predicted by prcomp object
#' @param x Character indicating PC to plot on x-axis
#' @param y Character indicating PC to plot on y-axis
#' @param ... Optional arguments are passed to aes_string in ggplot. Optional
#'   parameters have to match column names in metadata
ggplot_pca <- function(
  X,
  metadata,
  cex = 2,
  alpha = 1,
  show.legend = TRUE,
  plot.axis = TRUE,
  plot.rowname = FALSE,
  do.pca = TRUE,
  var_pc = NULL,
  newdata = NULL,
  x = "PC1", y = "PC2",
  ...
) {
  x_idx <- as.numeric(substring(x, 3))
  y_idx <- as.numeric(substring(y, 3))
  
  # PCA
  if (do.pca) {
    pca_obj <- prcomp(Matrix::t(X))
    Z <- data.frame(pca_obj$x)
    eigenvalues <- (pca_obj$sdev) ^ 2
    var_pc <- eigenvalues / sum(eigenvalues) # overwrites var_pc
    pc_names <- paste0("PC", seq_len(ncol(Z)))
    pc_labels <- sprintf("%s (%.2f%%)", pc_names, var_pc * 100)
  } else {
    Z <- data.frame(X)
    pc_names <- paste0("PC", seq_len(ncol(Z)))
    pc_labels <- sprintf("%s (%.2f%%)", pc_names, var_pc * 100)
  }
  # Projects newdata into PCA space
  if (!is.null(newdata)) {
    Z_new <- predict(pca_obj, newdata = Matrix::t(newdata))
    # Remove duplicate rows
    Z_new <- Z_new[!(rownames(Z_new) %in% rownames(Z)), ]
    Z <- rbind(Z, Z_new)
  }
  
  # Concat with metadata
  metadata_cols <- unlist(list(...))
  metadata1 <- metadata[rownames(Z), metadata_cols, drop = F]
  Z_metadata <- cbind(tibble::rownames_to_column(Z), metadata1)
  Z_metadata$rowname <- substring(Z_metadata$rowname, 1, 4)
      
  ax <- ggplot(
    Z_metadata,
    aes_string(x = x, y = y, label = "rowname", ...),
  ) +
    labs(x = pc_labels[x_idx], y = pc_labels[y_idx])

  if (plot.axis) {
    ax <- ax +
      geom_vline(xintercept = 0, color = "black") +
      geom_hline(yintercept = 0, color = "black")
  }
  # Plot text labels instead of points
  if (plot.rowname){
    ax <- ax +
      geom_point(
        color = "black", cex = cex,
        alpha = alpha, show.legend = show.legend
      ) +
      geom_text(cex = cex)
  } else {
    ax <- ax +
      geom_point(cex = cex, alpha = alpha, show.legend = show.legend)
  }

  ax
}


#' Plots top PCs using ggplot
#'
#' @import ggplot2
#' @importFrom tidyr gather
#' @export
#' @param X dataframe with features as rows and samples as columns
#' @param metadata dataframe of metadata with samples as rows
ggplot_top_pc <- function(
  X, metadata, x_axis, n = 8, nrow = 2,
  cex = 2, show.legend = FALSE, show.xaxis = TRUE,
  newdata = NULL, ...
) {
  # PCA
  pca_obj <- prcomp(t(X))
  if (ncol(pca_obj$x) < n) {
    n <- ncol(pca_obj$x)
    cat(sprintf("PCA transformed matrix only has %d dimensions!\n", n))
  }
  Z <- data.frame(pca_obj$x[, seq(n)])
  eigenvalues <- (pca_obj$sdev) ^2
  var_pc <- eigenvalues[seq(n)] / sum(eigenvalues)
  pc_labels <- sprintf("PC%d (%.2f%%)", seq(n), var_pc*100)
  names(pc_labels) <- paste0('PC', seq_len(length(pc_labels)))

  # Projects newdata into PCA space
  if (!is.null(newdata)) {
    Z_new <- predict(pca_obj, newdata = t(newdata))[, 1:3]
    # Remove duplicate rows
    Z_new <- Z_new[!(rownames(Z_new) %in% rownames(Z)), ]
    Z <- rbind(Z, Z_new)
  }

  # Concat with metadata
  plot_factors <- unlist(list(...))
  plot_factors <- unique(c(x_axis, plot_factors))
  metadata1 <- metadata[rownames(Z), plot_factors, drop = F]
  Z_metadata <- cbind(Z, metadata1)
  # Convert data to long format
  Z_long <- tidyr::gather(Z_metadata, key = "PC", value = "value", -plot_factors)
  
  ax <- ggplot(Z_long, aes_string(x = x_axis, y = "value", ...)) +
    facet_wrap(
      ~PC, scales = 'free_y', nrow = nrow,
      labeller = as_labeller(pc_labels),
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 1),
      cex = cex, alpha = 0.8,
      show.legend = show.legend
    )
  
  if (!show.xaxis) {
    ax <- ax +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  ax 
}


#' Plots UMAP plot using ggplot2
#'
#' @import ggplot2
#' @importFrom umap umap
#' @param obj UMAP object from umap function
#' @param ... optional arguments are passed to aes_string in ggplot
ggplot_umap <- function(
  X, metadata,
  cex = 2,
  plot_label = FALSE,
  return_data = FALSE, ...
) {
  obj <- umap(t(X))
  Z <- data.frame(obj$layout)
  colnames(Z) <- c("UMAP1", "UMAP2")
  
  metadata_cols <- unlist(list(...))
  metadata1 <- metadata[colnames(X), metadata_cols, drop = F]
  Z_metadata <- cbind(tibble::rownames_to_column(Z), metadata1)
  Z_metadata$rowname <- substring(Z_metadata$rowname, 1, 4)
  
  ax <- ggplot(
    Z_metadata,
    aes_string(x = 'UMAP1', y = 'UMAP2', label = "rowname", ...)
  ) 
  
  # Plot text labels instead of points
  if (plot_label) {
    ax <- ax +
      geom_point(color = "black", cex = cex, alpha = .2) +
      geom_text(cex = cex)
  } else {
    ax <- ax + geom_point(cex = cex)
  }

  if (return_data)
    return(list(plot = ax, X = Z_metadata))

  ax
}


#' Plots t-SNE plot using ggplot2
#'
#' @import ggplot2
#' @import Rtsne
#' @param obj UMAP object from umap function
#' @param ... optional arguments are passed to aes_string in ggplot
ggplot_tsne <- function(X, metadata, cex = 2, alpha = 1, ...) {
  obj <- Rtsne(normalize_input(t(X)))
  Z <- data.frame(obj$Y)
  colnames(Z) <- c("tSNE1", "tSNE2")
  metadata_cols <- unlist(list(...))
  metadata1 <- metadata[colnames(X), metadata_cols, drop = F]
  Z_metadata <- cbind(Z, metadata1)
  
  ggplot(Z_metadata, aes_string(x = "tSNE1", y = "tSNE2" , ...)) +
    geom_point(cex = cex, alpha = alpha)
}


#' Plots ROC curve
#'
#' @import ggplot2 pROC
#'
#' @param X dataframe containing scores / predictions and labels
#' @param response character containing name of column in X with labels
#' @param predictor character containing name of column/s in X of predictors
#' @param direction character belonging to c("auto", "<", ">"). control (0) < case (1)
#' no character vectors allowed to specific multiple different directions!
#' @param pauc.limits numeric vector of length 2 indicating limits
#' @param show.names logical indicating whether to include names in AUC legend
#' @param plot.names character containing names of each plot in same order as predictor.
#' Default NULL, predictor used for names
ggplot_roc <- function(
  X, response, predictor,
  direction = "auto",
  pauc.limits = FALSE,
  pauc.axis = c("specificity", "sensitivity"),
  pauc.correct = TRUE,
  plot.names = NULL,
  show.names = TRUE,
  return.auc = FALSE,
  lwd = 1
) {
  pauc.axis = match.arg(pauc.axis)
  
  auc_caption <- "AUC"
  if (!is.logical(pauc.limits)) {
    # If pAUC limits are provided
    pauc.limits <- sort(pauc.limits, decreasing = TRUE)
    auc_caption <- "pAUC"
  }
  
  if (length(predictor) == 1) {
    # Single ROC curve
    roc_objs <- pROC::roc(
      X[[response]],
      X[[predictor]],
      direction = direction,
      partial.auc = pauc.limits,
      partial.auc.focus = pauc.axis,
      partial.auc.correct = pauc.correct
    )
    aucs <- roc_objs$auc
    d <- data.frame(
      FPR = 1 - roc_objs$specificities,
      TPR = roc_objs$sensitivities
    )
    d <- d[nrow(d):1, ]
    d <- cbind(names = predictor, d)
  } else if (length(predictor) > 1) {
    # To ensure that direction can be subsetted
    if (length(direction) == 1) {
      direction <- rep(direction, length(predictor))
    } else if (length(predictor) != length(direction)) {
      stop("Length of direction not equals to that of predictor!")
    }
    
    roc_objs <- lapply(seq_along(predictor),
      function(i) {
        pROC::roc(
          X[[response]],
          X[[predictor[i]]],
          direction = direction[i],
          partial.auc = pauc.limits,
          partial.auc.focus = pauc.axis,
          partial.auc.correct = pauc.correct
        )
      }
    )
    aucs <- sapply(roc_objs, function(obj) obj$auc)
    list_d <- lapply(roc_objs, function(obj) data.frame(
      FPR = 1 - obj$specificities,
      TPR = obj$sensitivities
    ))
    list_d <- lapply(list_d, function(d) d[nrow(d):1, ])
    d <- do.call(rbind, list_d)
    n_rows <- sapply(list_d, nrow)
    predictors_col <- rep(predictor, n_rows)
    d <- cbind(names = predictors_col, d)
  } else {
    stop("arg predictor is of non-positive length.")
  }
        
  ## Plot labels
  if (show.names) {
    if (!is.null(plot.names)) {
      if (length(plot.names) != length(predictor))
        stop("length of plot names does not match number of predictors!")
      
      # plot.names in same order as predictor
      labels <- sprintf("%s (%s: %.3f)", plot.names, auc_caption, aucs)
      labels[is.na(aucs)] <- plot.names[is.na(aucs)]
    } else {
      labels <- sprintf("%s (%s: %.3f)", predictor, auc_caption, aucs)
      labels[is.na(aucs)] <- predictor[is.na(aucs)]
    }
  } else {
    labels <- sprintf("%s = %.3f", auc_caption, aucs)
    labels[is.na(aucs)] <- "pAUC = NA"
  }
  # Order according to lexical order of predictor
  labels <- labels[order(predictor)]

  ax_roc <- ggplot() +
    geom_segment(
      aes(x = 0, y = 0, xend = 1, yend = 1),
      inherit.aes = FALSE,
      lty = "dotted", lwd = lwd,
      colour = "black", alpha = .4
    ) +
    geom_line(
      data = d,  # to avoid mutiple plotting of geom_segment
      aes(x = FPR, y = TPR, col = names),  # OPTION lty = names
      direction = "hv", lwd = lwd
    ) +
    scale_color_discrete(
      name = element_blank(),
      label = labels
    ) +
#     scale_linetype_discrete(
#       name = element_blank(),
#       label = labels
#     ) +
    theme_bw() +
    labs(x = "FPR", y = "TPR") +
    theme(
      legend.position = c(.95, .05),
      legend.justification = c("right", "bottom"),
      legend.background = element_rect(fill = NA)
    )
  
  if (is.logical(pauc.limits)) {
    # if no pauc.limits is provided
    ax_roc <- ax_roc +
      coord_cartesian(xlim = c(0, 1)) +
      coord_cartesian(ylim = c(0, 1))
  } else if (pauc.axis == "specificity") {
    fpr_limits <- 1 - pauc.limits
    d_rect <- data.frame(
      xmin = fpr_limits[1], xmax = fpr_limits[2],
      ymin = 0, ymax = 1
    )
    ax_roc <- ax_roc +
#       coord_cartesian(xlim = fpr_limits) +
      geom_rect(
        data = d_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "blue", alpha = 0.2
      )
  } else if (pauc.axis == "sensitivity") {
    d_rect <- data.frame(
      xmin = 0, xmax = 1,
      ymin = pauc.limits[2], ymax = pauc.limits[1]
    )
    ax_roc <- ax_roc +
#       coord_cartesian(ylim = pauc.limits) +
      geom_rect(
        data = d_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "blue", alpha = 0.2
      )
  }
  
  if (return.auc) {
    return(list(plot = ax_roc, auc = aucs))
  }
  ax_roc
}


# Plots ROC and calculates AUC in a primitive fashion (i.e. ROC is step function)
# Does not resolve ties in the score
# Assumption: Lower score will be labelled preferentially as 1, ROC is step function
# Assumption that score vec and label vecs are corresponding
plot_roc <- function(score_list, label_vec,
                    name_vec = NULL, plot_title = NULL,
                    is_bigger_better_vec = rep(F, length(score_list))) {
  # Function to plot a single ROC curve and calculate AUC
  plotSingleROC <- function(score_vec, is_bigger_better, color) {
    # Sort label vector according to score vector in ascending order
    sort_label <- label_vec[order(score_vec, decreasing = is_bigger_better)]
    # Dataframe of TPR and FPR
    df <- data.frame(TPR = cumsum(sort_label)/sum(sort_label),
                     FPR = cumsum(!sort_label)/sum(!sort_label))
    # Insert 0, 0 at the first row
    roc_df <- rbind(c(0,0), df)
    # Calculates change in FPR at each step
    dFPR <- c(0, diff(roc_df$FPR))
    # Sum of individual rectangle steps
    AUC <- sum(roc_df$TPR * dFPR)
    # Plot single ROC curve
    lines(roc_df$FPR, roc_df$TPR,
          col = color, lwd = 3)
    return(AUC)
  }
  # Initialise plot
  colour_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                      "#66A61E", "#E6AB02", "#A6761D", "#666666")
  color_index <- colour_palette[1:length(score_list)]
  plot(NULL, main = plot_title,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0,1), ylim = c(0,1),
       xaxs = "i", yaxs = "i",
       cex.main = 1.8, cex.lab = 1.7, cex.axis = 1.6)
  abline(0, 1, lty = 5, lwd = 2)
  auc_vec <- mapply(plotSingleROC, score_list, is_bigger_better_vec,
                    color_index, SIMPLIFY = T)
  # If name_vec is not NULL display legend
  if (!is.null(name_vec)) {
    format_text <- function(name, auc) sprintf("%s (%.3f)", name, auc)
    legend_text <- mapply(format_text, name_vec, auc_vec)
    legend("bottomright", inset = 0.03, lty = 1, lwd = 2,
           cex = 1.2, bty = "o", text.width = 0.35,
           legend = legend_text, col = color_index)
  }
  return(auc_vec)
}


#' Performs PCA and plots vectors. Assume: 3 normal samples
#' @param X dataframe with features x samples and containing
#' D0, D8 and 3 Normal patients (ordered correctly)
plot_vectors <- function(
  X, metadata, pca = T, cex = 3, main = NULL
) {
  # PCA
  if (pca) {
    pca_obj <- prcomp(t(X))
    pca_df <- data.frame(pca_obj$x[,1:2])
    n <- nrow(pca_df)
    d0_pca <- pca_df[1:((n-3) / 2),]
    d8_pca <- pca_df[((n-1) / 2):(n-3),]
    
    stopifnot(all(substring(rownames(d0_pca), 1, 4) == substring(rownames(d8_pca), 1, 4)))
    
    subtype_pca <- cbind(d0_pca, d8_pca)
    colnames(subtype_pca) <- c("start_x", "start_y", "end_x", "end_y")
    norm_pca <- pca_df[(n-2):n,]
    
    # Obtaining batch and class annotations
    label <- as.factor(metadata[rownames(d0_pca), "label"])
    batch <- as.factor(metadata[rownames(d0_pca), "batch_info"])
    
    # Axis labels
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:4]/sum(eigenvalues)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:4, var_pc*100)
  } else {
    print("No PCA performed!")
    pca_df <- data.frame(X)
  }
  
  scatter_pca <- ggplot(data = subtype_pca) +
    geom_point(
      aes(x = start_x, y = start_y, fill = batch), 
      shape = 21, size = cex, show.legend = T
    ) +
    geom_point(
      aes(x = end_x, y = end_y, fill = batch),
      shape = 22, size = cex, show.legend = F
    ) +
    geom_segment(
      aes(x = start_x, y = start_y, xend = end_x, yend = end_y, colour = label),
      arrow = arrow(length = unit(0.3, "cm")), alpha = 0.5
    ) +
#     geom_text(aes(x = PC1_A, y = PC2_A,
#                   label = rownames(subtype_pca)),
#               position = position_nudge(x = 4, y = 2), size = 2.5) +
    geom_point(
      data = norm_pca,
      aes(x = PC1, y = PC2),
      size = cex, shape = 17
    ) +
    scale_color_manual(values = COL_LABEL)
  
  if (pca) {
    scatter_pca <- scatter_pca +
      xlab(pc_labels[1]) + ylab(pc_labels[2])
  }
  
  if (!is.null(main)) {
    scatter_pca <- scatter_pca + labs(title = main)
  }
  
  return(scatter_pca)
}


#' Plots boxplot of features
#' Provides p-values from wilcoxon rank-sum test
#' If different group and color is provided X_y is flattened accordingly
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @export
plot_boxplot <- function(
  X_y,
  features = c("erm1_ratio2", "l2norm_ratio2", "angle_LD0_LD8_ratio2", "log_mrd_d33"),
  group = 'label',
  fill = 'label',
  p.value = TRUE,
  show.legend = FALSE,
  ...
) {
  aes_extra <- list(...)
  pch_treatment <- 21:24
  names(pch_treatment) <- c('SR', 'IR', 'HR1', 'HR2')
  
  # Assert that all features are present in data provided
  feature_idx <- sapply(features, match, colnames(X_y))
  stopifnot(!is.na(sum(feature_idx)))
  
  # WARNING: Manual inputted values 
  feature_order <- c(features, "p_d8", "p_d33")
  feature_labels <- c(
    "'ERM Ratio'", "'ARM Ratio'", "phi",
    "log[10](MRD)", "paste('P(Remission|', bold(x['D8']), ', s)')",
    "paste('P(Remission|', bold(x['D33']), ', s)')"
  )
 
  X_y <- tibble::rownames_to_column(X_y)
  features <- colnames(X_y)
  # All features aside from those in feature_order are gathered
  patient_info <- unique(setdiff(features, feature_order))
  long_X_y <- gather(X_y, key = "feature", value = "value", -patient_info)
  
  # Reorder levels and label features
  long_X_y$feature <- factor(
    long_X_y$feature,
    levels = feature_order,
    labels = feature_labels
  )
  
  ax_jitter <- ggplot(long_X_y) +
    facet_wrap(
      ~feature,
      nrow = 1, scales = "free",
      labeller = label_parsed
    ) +
    geom_boxplot(
      aes_string(x = group, y = 'value'),
      col = "black", alpha = 0,
      show.legend = show.legend
    ) +
    # geom_text(
    #   aes(col = label),
    #   position = position_jitterdodge(jitter.width = 1, seed = 1),
    #   cex = 2.5, show.legend = FALSE
    # ) + 
    scale_shape_manual(values = pch_treatment) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 15),
      axis.title.y = element_blank()
      # legend.position = "none"
    )

  ### Plot options ###
  if (fill == 'label') {
    ax_jitter <- ax_jitter + scale_fill_manual(values = COL_LABEL)
  }

  # if user specificies pch
  if ('pch' %in% names(aes_extra)) {
    ax_jitter <- ax_jitter +
      geom_point(
        aes_string(x = group, y = 'value', fill = fill, ...),
        position = position_jitterdodge(jitter.width = 1, seed = 1),
        cex = 2.5, alpha = 1, col = "black", show.legend = show.legend
      )
  } else {
    ax_jitter <- ax_jitter +
      geom_point(
        aes_string(x = group, y = 'value', fill = fill, ...),
        position = position_jitterdodge(jitter.width = 1, seed = 1),
        pch = 21, cex = 2.5, alpha = 1, col = "black",
        show.legend = show.legend
      )
  }

  # Jitter plot: p-value label
  has_classes <- length(table(X_y[[group]])) > 1
  # Both group sizes must be > 1
  has_samples <- min(table(X_y[[group]])) > 1
  if (p.value && has_classes && has_samples) {
    stat_label <- sapply(feature_order, function(idx) {
      x1_x2 <- split(X_y[, idx], X_y[[group]])
      tryCatch(
        {
          # Wilcoxon test p-value
          x1 <- x1_x2[[1]]
          x2 <- x1_x2[[2]]
          htest <- wilcox.test(x1, x2, exact = T)
          mean_diff <- mean(x1) - mean(x2)
          # Cohen's d assuming unequal var. Hence sd = sqrt((s_1^2 + s_2^2)/2)
          sd_unequal <- mean(c(var(x1), var(x2))) ^ .5
          cohen_d <- mean_diff / sd_unequal
          sprintf("p = %.3f\nd = %.2f", htest$p.value, cohen_d)
        },
        error = function(err) {
          print(err)
          return(err)
        }
      )
    })
   
    x_coord <- rep(1.55, 6)
    x_coord[4] <- 0.55

    compute_ycoord <- function(x) 0.92 * (max(x) - min(x)) + min(x)
    y_coord <- c(
      compute_ycoord(X_y[feature_order[1]]),
      compute_ycoord(X_y[feature_order[2]]),
      compute_ycoord(X_y[feature_order[3]]),        
      compute_ycoord(X_y[feature_order[4]]),
      compute_ycoord(X_y[feature_order[5]]),
      compute_ycoord(X_y[feature_order[6]])
    )
    
    ann_text <- data.frame(
      feature = factor(feature_order, levels = feature_order, labels = feature_labels),
      x_coord = x_coord,
      y_coord = y_coord, 
      stat_label = stat_label
    )
    ann_text <- na.omit(ann_text)
    
    ax_jitter <- ax_jitter +
      geom_text(
        data = ann_text,
        aes(x = x_coord, y = y_coord, label = stat_label),
        size = 3, hjust = 0 # colour = "black"
      )
  }
  
  ax_jitter
}


#' Plots boxplots according to colour and shape provided
#' No statistical tests conducted between groups.
#' Plots both proba (D8) and (D33)
#'
#' @param colour_by string indicating name feature to group colour by
#' @import ggplot2
#' @importFrom tidyr gather
#' @export
plot_boxplots_v2 <- function(X_y, colour_by, shape_by) {
  FEAT_ORDER <- c(
    "erm1_ratio2", "l2norm_ratio2", "angle_d0d8_d0normal",
    "log_mrd_d33", "p_d8", "p_d33"
  )
  FEAT_LABS <- c(
    "'ERM Ratio'", "'ARM Ratio'", "theta",
    "log[10](MRD)", "paste('P(Remission|', bold(x[D8]), ', s)')",
    "paste('P(Remission|', bold(x[D33]), ', s)')"
  )
 
  long_X_y <- gather(
    X_y, key = "feature", value = "value",
    -c(colour_by, shape_by)
  )
  # Reorder levels and label features
  long_X_y$feature <- factor(
    long_X_y$feature,
    levels = FEAT_ORDER,
    labels = FEAT_LABS
  )
  
  ax_jitter <- ggplot(
    long_X_y,
    aes_string(x = colour_by, y = "value")
  ) +
    # geom_boxplot(
    #   aes_string(group = colour_by),
    #   col = "black", alpha = 0,
    #   show.legend = F
    # ) +
    geom_point(
      aes_string(fill = colour_by, shape = shape_by),
      position = position_jitterdodge(jitter.width = 1),
      cex = 2.5, alpha = 1, show.legend = F
    ) +
    scale_shape_manual(values = 21:22) +
    facet_wrap(
      ~feature,
      nrow = 1, scales = "free",
      labeller = label_parsed
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 15),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  ax_jitter
}


#' Plots boxplot of features
#' Provides p-values from wilcoxon rank-sum test
#' If different group and color is provided X_y is flattened accordingly
#'
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @export
plot_boxplot1 <- function(
  X_y,
  features = c("erm1_ratio2", "l2norm_ratio2", "angle_LD0_LD8_ratio2", "log_mrd_d33"),
  group = 'label',
  fill = 'label',
  p.value = TRUE,
  show.legend = FALSE,
  ...
) {
  aes_extra <- list(...)
  pch_treatment <- 21:24
  names(pch_treatment) <- c('SR', 'IR', 'HR1', 'HR2')
  
  # Assert that all features are present in data provided
  feature_idx <- sapply(features, match, colnames(X_y))
  stopifnot(!is.na(sum(feature_idx)))
  
  # WARNING: Manual inputted values 
  feature_order <- c(features, "p_d8", "p_d33")
  feature_labels <- c(
    "'ERM Ratio'", "'ARM Ratio'", "phi",
    "log[10](MRD)", "paste('P(Remission|', bold(x['D8']), ', s)')",
    "paste('P(Remission|', bold(x['D33']), ', s)')"
  )
 
  X_y <- tibble::rownames_to_column(X_y)
  features <- colnames(X_y)
  # All features aside from those in feature_order are gathered
  patient_info <- unique(setdiff(features, feature_order))
  print(patient_info)
  long_X_y <- gather(X_y, key = "feature", value = "value", -patient_info)
  print(long_X_y) 

  # Reorder levels and label features
  long_X_y$feature <- factor(
    long_X_y$feature,
    levels = feature_order,
    labels = feature_labels
  )
  
  ax_jitter <- ggplot(
    long_X_y,
    aes_string(x = group, y = 'value', label = 'rowname', fill = fill)
  ) +
    # geom_boxplot(
    #   aes_string(group = group),
    #   col = "black", alpha = 0,
    #   show.legend = show.legend
    # ) +
    # geom_text(
    #   aes(col = label),
    #   position = position_jitterdodge(jitter.width = 1, seed = 1),
    #   cex = 2.5, show.legend = FALSE
    # ) + 
    scale_shape_manual(values = pch_treatment) +
    facet_wrap(
      ~feature,
      nrow = 1, scales = "free",
      labeller = label_parsed
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 15),
      axis.title.y = element_blank()
      # legend.position = "none"
    )

  ### Plot options ###
  if (fill == 'label') {
    ax_jitter <- ax_jitter + scale_fill_manual(values = COL_LABEL)
  }

  if ('pch' %in% names(aes_extra)) {
    ax_jitter <- ax_jitter +
      geom_point(
        aes_string(...),
        position = position_jitterdodge(jitter.width = 1, seed = 1),
        cex = 2.5, alpha = 1, col = "black", show.legend = show.legend
      )
  } else {
    ax_jitter <- ax_jitter +
      geom_point(
        aes_string(...),
        position = position_jitterdodge(jitter.width = 1, seed = 1),
        pch = 21, cex = 2.5, alpha = 1, col = "black",
        show.legend = show.legend
      )
  }

  # Jitter plot: p-value label
  has_classes <- length(table(X_y[[group]])) > 1
  # Both group sizes must be > 1
  has_samples <- min(table(X_y[[group]])) > 1
  if (p.value && has_classes && has_samples) {
    
    stat_label <- sapply(feature_order, function(idx) {
      x1_x2 <- split(X_y[, idx], X_y[[group]])
      tryCatch(
        {
          # Wilcoxon test p-value
          x1 <- x1_x2[[1]]
          x2 <- x1_x2[[2]]
          htest <- wilcox.test(x1, x2, exact = T)
          mean_diff <- mean(x1) - mean(x2)
          # Cohen's d assuming unequal var. Hence sd = sqrt((s_1^2 + s_2^2)/2)
          sd_unequal <- mean(c(var(x1), var(x2))) ^ .5
          cohen_d <- mean_diff / sd_unequal
          sprintf("p = %.3f\nd = %.2f", htest$p.value, cohen_d)
        },
        error = function(err) {
          print(err)
          return(err)
        }
      )
    })
   
    x_coord <- rep(1.55, 6)
    x_coord[4] <- 0.55

    compute_ycoord <- function(x) 0.92 * (max(x) - min(x)) + min(x)
    y_coord <- c(
      compute_ycoord(X_y[feature_order[1]]),
      compute_ycoord(X_y[feature_order[2]]),
      compute_ycoord(X_y[feature_order[3]]),        
      compute_ycoord(X_y[feature_order[4]]),
      compute_ycoord(X_y[feature_order[5]]),
      compute_ycoord(X_y[feature_order[6]])
    )
    
    ann_text <- data.frame(
      feature = factor(feature_order, levels = feature_order, labels = feature_labels),
      x_coord = x_coord,
      y_coord = y_coord, 
      stat_label = stat_label
    )
    ann_text <- na.omit(ann_text)

    ax_jitter <- ax_jitter +
      geom_text(
        data = ann_text,
        aes(x = x_coord, y = y_coord, label = stat_label),
        size = 3, hjust = 0 # colour = "black"
      )
  }
  
  ax_jitter
}


# 3D PCA plot
rglplot_scatter <- function(
  df, colour = "blue", pch = 21, pc_labels = NULL,
  ratio_list = list(2,1,1)
) {
  # RGL plot parameters
  rgl.open()
  rgl.bg(color="white")
  rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  par3d(windowRect = c(50, 20, 500, 500))
  pch3d(df[,1], df[,2], df[,3], bg = colour,
        pch = pch, cex = 0.5, lwd = 1.5)
  box3d(col = "black")
  # title3d(xlab = pc_labels[1], ylab = pc_labels[2],
  #         zlab = pc_labels[3], col = "black")
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}


# 3D PCA plot
plot3d_pca <- function(
  df, colour, pch, pc_labels = NULL,
  cex = 0.5, lwd = 0.5,
  width = 600, height = 400,
  ratio_list = list(2, 1, 1)
) {
  if (is.null(pc_labels)) {
    print("PCA performed!")
    pca_obj <- prcomp(t(df)) 
    pca_df <- as.data.frame(pca_obj$x[,1:3])
    eigenvalues <- (pca_obj$sdev)^2
    var_pc <- eigenvalues[1:3]/sum(eigenvalues)
    pc_labels <- sprintf("PC%d (%.2f%%)", 1:3, var_pc*100)
  } else {
    print("No PCA performed!")
    pca_df <- as.data.frame(df)
  }
  
  # RGL plot parameters
  if (rgl.cur() == 0) {
    rgl.open()
    rgl.bg(color = "white")
  }
  rgl.clear(type = "all")  
  # rgl.viewpoint(zoom = 0.8)
  # rgl.viewpoint(theta = 110, phi = 5, zoom = 0.8)
  
  # figure size
  window <- c(0, 0, width, height) + 50
  par3d(windowRect = window) 
  
  with(pca_df, pch3d(
    PC1, PC2, PC3,
    bg = colour, pch = pch,
    cex = cex, lwd = 0.5
  ))
  box3d(col = "black")
  title3d(
    xlab = pc_labels[1], ylab = pc_labels[2],
    zlab = pc_labels[3], col = "black"
  )
  
  # Plot aspect ratios of axis according to variance
  do.call(aspect3d, ratio_list)
}


# Plot PCA before selecting features
# Batch information of all the timepoints
plot3d_yeoh <- function(X, metadata, ...) {
  # colour
  batch <- metadata[colnames(X), "batch_info"]
  batch_levels <- sort(unique(batch))
  batch_palette <- ggplot_palette(length(batch_levels))
  names(batch_palette) <- batch_levels
  sample_colour <- batch_palette[as.character(batch)]
  # shape
  class <- metadata[colnames(X), "class_info"]
  levels(class) <- 21:23
  sample_shape <- as.numeric(as.character(class))
  
  plot3d_pca(X, sample_colour, sample_shape, ...)
}


plotHeatmapSubtype <- function(X, metadata_df) {
  subtype_factor <- metadata_df[colnames(X), "subtype"]
  set3_pal <- brewer.pal(9, "Set3")
  subtype_col <- set3_pal[subtype_factor]
  
  par(mar = c(1,1,1,1))
  heatmap(data.matrix(X),
          col = brewer.pal(9, "Blues"),
          ColSideColors = subtype_col,
          scale = "none",
          labRow = NA, labCol = NA)
  
  legend(x = -.04, y = 1350, legend = levels(subtype_factor),
         col = set3_pal[factor(levels(subtype_factor))],
         pch = 15, cex = .7)
  heatmap_subtype <- recordPlot()
  par(mar = c(5.1, 4.1, 4.1, 2.1)) # Reset to defaults
  return(heatmap_subtype)
}

#' @import dplyr
plot_mean <- function(df, batch_vec1) {
  # Melt dataframe
  melt_df <- gather(df, key = "ID", value = "value")
  print(head(melt_df))
  # Trimmed mean probe intensities for each chip
  mean_tibble <- melt_df %>%
    group_by(ID) %>%
    summarise(mean = mean(value))
  mean_batch_tibble <- cbind(mean_tibble,
                             batch_vec1 = batch_vec1[mean_tibble$ID])
  
  mean_scatter <- ggplot(mean_batch_tibble, aes(x = ID, y = mean)) +
    geom_point(aes(col = factor(batch_vec1)),
               show.legend = F, size = 3) +
    facet_wrap(factor(batch_vec1), scales = "free_x") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  return(mean_scatter)
}


# Assumes that dataframe has been log-transformed
plotExplore <- function(df1, metadata_df) {
  # Obtaining batch and class annotations
  batch_factor <- as.factor(metadata_df[colnames(df1),"batch_info"])
  class_factor <- metadata_df[colnames(df1),"class_info"]
  
  # Melt dataframe
  melt_df <- gather(df1, key = "ID", value = "value") 
  melt_df$batch <- as.factor(metadata_df[melt_df$ID,"batch_info"])
  melt_df$class <- metadata_df[melt_df$ID,"class_info"]
  
  # Plot means
  mean_tibble <- melt_df %>% group_by(ID) %>% summarise(mean = mean(value))
  mean_scatter <- ggplot(mean_tibble, aes(x = ID, y = mean,
                                          col = batch_factor,
                                          pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot boxplot
  boxplot <- ggplot(melt_df, aes(x = ID, y = value, col = ID)) + 
    geom_boxplot(show.legend = F) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Plot density curve
  print(head(melt_df))
  print(tail(melt_df))
  pdf <- ggplot(melt_df, aes(x = value, group = ID, col = batch_info)) +
    geom_density(show.legend = F, alpha = 0.3) +
    facet_wrap(~class_info) +
    scale_color_viridis_d()
  
  # Plot PCA
  # Filters out rows with all zero values
  nonzero_logvec <- rowSums(df1) != 0 & apply(df1, 1, var) != 0
  pca_obj <- prcomp(t(df1[nonzero_logvec,]))
  pca_df <- data.frame(pca_obj$x[,1:4])
  eigenvalues <- (pca_obj$sdev)^2
  var_pc <- eigenvalues[1:5]/sum(eigenvalues)
  pc_labels <- sprintf("PC%d (%.2f%%)", 1:5, var_pc*100)
  
  pc1_pc2 <- ggplot(pca_df, aes(x = PC1, y = PC2, col = batch_factor,
                                pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[2])
  pc1_pc3 <- ggplot(pca_df, aes(x = PC1, y = PC3, col = batch_factor,
                                pch = class_factor)) +
    geom_point(size = 3, show.legend = F) +
    labs(x = pc_labels[1], y = pc_labels[3])
  
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc1_pc3)
  multiplot <- plot_grid(mean_scatter, pdf, pca, nrow = 3)
  return(multiplot)
}


# Plot venn diagram
jacc_coeff <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                                  category = c("D1", "D2"),
                                  cex = 3, fontfamily = "sans",
                                  cat.cex = 3, cat.fontfamily = "sans",
                                  margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  print(unname(venn_area[3]/union))
  return(venn_plot)
}


#' Get column names in same order as pheatmap plot
#' @param obj pheatmap object
get_colnames <- function(obj) {
  obj$tree_col$labels[obj$tree_col$order]
}


#' Get row names in same order as pheatmap plot
#' @param obj pheatmap object
get_rownames <- function(obj) {
  obj$tree_row$labels[obj$tree_row$order]
}


#' Generate default ggplot colours
ggplot_palette <- function(n) {
  hues = seq(15, 375, length = n + 1)
  return(hcl(h = hues, c = 100, l = 65)[1:n])
}

#' @param D dataframe containing columns: c(score, label, treatment)
#' @param treatment character indicating treatment or treatment_current
#' to use for scoring and as fill color
plot_barchart <- function(
  D, timepoint,
  treatment = c("treatment", "treatment_current"),
  ylim = c(0, 35),
  show.legend = F
) {
  treatment <- match.arg(treatment)
  treatment_caption <- ifelse(
    treatment == "treatment", "Final treatment", "Current treatment"
  )
  # Fill colours
  colors <- c("dodgerblue", "chartreuse3", "tomato")
  risk_levels <- c("SR", "IR", "HR")
  names(colors) <- risk_levels
  
  # Awarding theoretical scores
  D1 <- D[c("label", treatment)] # final or current treatment
  colnames(D1)[2] <- "treatment"
  theoretical_scoring_table <- data.frame(
    label = factor(
      rep(c("Remission", "Relapse"), each = 3),
      levels = c("Remission", "Relapse")
    ),
    treatment = factor(rep(risk_levels, 2), levels = risk_levels),
    score = c("1", "1", "1", "0", "0", "1")
  )
  D2 <- merge(D1, theoretical_scoring_table, by = c("label", "treatment"))
  score_labs <- c("Score: 0", "Score: 1")
  names(score_labs) <- c("0", "1") 

  ggplot(D) +
    facet_grid(
      label ~ score,
      labeller = labeller(score = score_labs)
    ) +
    geom_bar(
      aes_string(x = "prediction", fill = treatment), # final or current treatment
      show.legend = show.legend
    ) +
    geom_bar(
      data = D2, aes(x = treatment),
      alpha = 0, color = "black", show.legend = show.legend
    ) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = sprintf(
        "%s (Total: %d/%d)",
        timepoint, sum(D$score), nrow(D)
      ),
      x = "Predicted risk level", y = "Count",
      fill = treatment_caption 
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 15)
    ) +
    ylim(ylim)
}


#' Plots cumulative sum of squares and HVP across features
#'
#' @param res result object from HVP 
#' @param m numeric indicating the top number of features according to sum of
#'   squares to plot. Defaults to plotting all features.
#'
#' @returns A ggplot object
#'
#' @import ggplot2 cowplot
#'
#' @keywords internal
#' @noRd
#'
plot_hvp <- function(res, m = NULL) {
  xlab <- "Feature index"
  # sort by sum of squares total in decreasing order
  ss <- res$sum.squares[rev(order(res$sum.squares[, "ss_total"])), ]
  HVP <- res$HVP

  if (!is.null(m))
    ss <- ss[seq_len(m), , drop = FALSE]

  ax_ssb <- ggplot(ss) +
    geom_line(
      aes(x = seq_len(nrow(ss)), y = cumsum(ss_batch)),
      color = "blue"
    ) +
    labs(x = xlab, y = "Cumulative sum of squares batch")

  ax_sst <- ggplot(ss) +
    geom_line(
      aes(x = seq_len(nrow(ss)), y = cumsum(ss_total)),
      color = "chartreuse3"
    ) +
    labs(x = xlab, y = "Cumulative sum of squares total")

  ax_hvp <- ggplot(ss) +
    geom_line(
      aes(x = seq_len(nrow(ss)), y = cumsum(ss_batch) / cumsum(ss_total)),
      col = "brown1"
    ) +
    geom_hline(yintercept = HVP, col = "darkgray", linetype = "dashed") +
    labs(x = xlab, y = "Cumulative HVP") +
    ylim(0, 1)

  plot_grid(ax_ssb, ax_sst, ax_hvp, nrow = 1)
}
