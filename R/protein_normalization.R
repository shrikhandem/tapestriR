
#' Log normalization factor per feature
#'
#' @param x Sample name
#'
#' @return normalization factor
#' @export
#' @examples
#' \dontrun{
#' norm_factor <- log_colMeans(data)
#' }
log_colMeans <- function(x){
  if(is.vector(x)) x <- matrix(x, nrow=1)
  return(exp(colMeans(log(x))))
}


#' Log normalization by cell
#'
#' @param counts counts data
#'
#' @return normalized counts
#' @export
#' @examples
#' \dontrun{
#' counts.norm <- log_norm_by_cell(counts)
#' }
log_norm_by_cell <- function(counts) {
  counts = counts + 1
  norm_counts <- (log(counts /colMeans(counts)))
  return(norm_counts)
}


#' Log normalization by feature
#'
#' @param counts counts data
#'
#' @return normalized counts
#' @export
#' @examples
#' \dontrun{
#' counts.norm <- log_norm_by_feature(counts)
#' }
log_norm_by_feature <- function(counts) {
  norm_counts <- counts + 1

  norm_counts <- t(apply(
    norm_counts, 1,
    norm <- function(x) {
      return(log(x / mean(as.numeric(x), na.rm = FALSE)))
    }
  ))

  return(norm_counts)
}


#' CLR normalization by cell
#'
#' @param counts counts data
#'
#' @return normalized counts
#' @export
#' @importFrom  SciViews ln
#' @importFrom  EnvStats geoMean
#' @examples
#' \dontrun{
#' counts.norm <- clr_by_feature(counts)
#' }
clr_by_cell <- function(counts) {
  counts = counts + 1
  clr_counts <- (SciViews::ln(counts /log_colMeans(counts)))
  return(clr_counts)
}


#' CLR normalization by feature
#'
#' @param counts counts data
#'
#' @return normalized counts
#' @export
#' @importFrom  SciViews ln
#' @importFrom  EnvStats geoMean
#' @examples
#' \dontrun{
#' counts.norm <- clr_by_feature(counts)
#' }
clr_by_feature <- function(counts) {
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

