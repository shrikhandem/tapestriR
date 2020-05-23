
# divide by the columns
div_col <- function(df, z, na.rm = FALSE) {
  df <- as.matrix(df)
  divnorm <- c()
  for (col in colnames(df)) {
    divnorm <- cbind(divnorm, as.numeric(round(df[, col] / z[col], 3)))
  }
  colnames(divnorm) <- colnames(df)
  rownames(divnorm) <- rownames(df)
  return(divnorm)
}

# Log normalization factor per feature
#
# @param x Sample name
# @return normalization factor
# @examples
# \dontrun{
# norm_factor <- log_colMeans(data)
# }
log_colMeans <- function(x){
  if(is.vector(x)) x <- matrix(x, nrow=1)
  return(exp(colMeans(log(x))))
}


# Log normalization by cell
#
# @param counts counts data
# @return normalized counts
# @examples
# \dontrun{
# counts.norm <- log_norm_by_cell(counts)
# }
log_norm_by_cell <- function(counts, ...) {
  counts = counts + 1
  norm_counts <- (log(counts /colMeans(counts)))
  return(norm_counts)
}


# Log normalization by feature
#
# @param counts counts data
# @return normalized counts
# @examples
# \dontrun{
# counts.norm <- log_norm_by_feature(counts)
# }
log_norm_by_feature <- function(counts, ...) {
  norm_counts <- counts + 1

  norm_counts <- t(apply(
    norm_counts, 1,
    norm <- function(x) {
      return(log(x / mean(as.numeric(x), na.rm = FALSE)))
    }
  ))

  return(norm_counts)
}


# CLR normalization by cell
#
# @param counts counts data
# @return normalized counts
# @examples
# \dontrun{
# counts.norm <- clr_by_cell(counts)
# }
clr_by_cell <- function(counts, ...) {
  counts = counts + 1
  clr_counts <- (SciViews::ln(counts /log_colMeans(counts)))
  return(clr_counts)
}


# CLR normalization by feature
#
# @param counts counts data
# @return normalized counts
# @examples
# \dontrun{
# counts.norm <- clr_by_feature(counts)
# }
clr_by_feature <- function(counts, ...) {
  counts = counts + 1
  clr_counts_wCells <- counts
  clr_counts_wCells <- t(apply(
    clr_counts_wCells, 1,
    norm <- function(x) {
      return(SciViews::ln(x / EnvStats::geoMean(as.numeric(unlist(c(x))), na.rm = FALSE)))
    }
  ))
  return(clr_counts_wCells)
}


#' Normalize raw counts counts
#'
#' @param analyte Analyte object to be normalized (Supported analytes : prot/dna)
#' @param group_by normalization to be performed across rows or columns (Supported : cell/feature)
#' @param method Normalization method (Supported : mean/median/clr/sum)
#' @return Analyte Analyte object with normalized data
#' @export
#' @examples
#' \dontrun{
#' dna_analyte <- normalize_data(dna_analyte, group_by="cell", method = "median")
#' prot_analyte <- normalize_data(prot_analyte, group_by="cell", method = "clr")
#' }
normalize_data <- function(analyte, group_by, method) {
  err.msg1 <- "Normalization method %s grouped by %s not supported for %s. Please check documentation for usage."
  err.msg2 <- "Group by %s not supported for %s. Please check documentation for usage."
  methods_supported <- c("median", "mean", "clr", "sum")
  assaytype <- analyte$type
  counts <- analyte$counts
  advanced <- FALSE
  assay_supported <- c("prot", "dna")
  if (!(assaytype %in% assay_supported)) {
    stop("Please use valid assay type (dna/prot)")
  }
  if (!(method %in% methods_supported)) {
    stop("Please provide valid method (median/mean/clr/sum)")
  }

  func <- NULL

  # Normalization methods from prot analyte
  if (assaytype == "prot") {
    if (method == "clr" && group_by == "feature") {
      func <- clr_by_feature
    }
    else if (method == "clr" && group_by == "cell") {
      func <- clr_by_cell
    } 
    else if (method == "mean" && group_by == "feature") {
      func <- log_norm_by_feature
    }
    else if (method == "mean" && group_by == "cell") {
      func <- log_norm_by_cell
    }
    else if (method == "sum" && group_by == "cell") {
      func <- normalize_barcodes
    } else {
      stop(sprintf(err.msg1, method, group_by, assaytype))
    } 
  }
  else if (assaytype == "dna") {
    if (group_by == "cell") {
      func <- normalize_barcodes
      if (method == "median") {
        advanced <- TRUE
      } else if (method == "sum") {
        advanced <- FALSE
      } else {
        stop(sprintf(err.msg1, method, group_by, assaytype))
      }
    } else {
      stop(sprintf(err.msg2, group_by, assaytype))
    }
  } else {
    stop(sprintf("Assay type %s not supported currently", assaytype))
  }

  if (is.null(func)) {
    stop("Error in running normalization. Please check documentation for usage.")
  }

  norm_data <- do.call(func, list(counts=counts, advanced=advanced))
  analyte$data <- data.frame(norm_data)
  analyte$normalized <- TRUE
  analyte$features <- colnames(norm_data)
  analyte$counts <- analyte$counts[, colnames(norm_data)]
  return(analyte)
}
