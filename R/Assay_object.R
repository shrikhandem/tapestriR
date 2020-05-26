
setOldClass(c("tbl_df", "tbl", "data.frame"))

#' Assay class
#'
#' Store assay specific details for given sample
#'
#' @slot assay_name 
#' @slot cell_annotations 
#' @slot feature_annotations 
#' @slot data_layers 
#' @slot analysis_layers 
#'
#' @name Assay-class
#' @rdname Assay-class
#' @import dplyr
#' @exportClass Assay
Assay <- setClass(
  Class = "Assay",
  slots = c(
    assay_name = "character",
    feature_annotations = "tbl_df",
    cell_annotations = "tbl_df",
    data_layers = "list",
    analysis_layers = "list"
  )
  #,contains = c(class(tibble()))
)
setClassUnion(name="AssayorNULL", c("Assay", "NULL"))

#' @export
#' @method dim Assay
#'
dim.Assay <- function(x) {
  c(nrow(x@cell_annotations), nrow(x@feature_annotations))
}

#' @export
#' @method nrow Assay
#'
nrow.Assay <- function(x) {
  nrow(x@cell_annotations)
}


#' @export
#'
"$.Assay" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#' @export
#'
"$<-.Assay" <- function(x, i, ..., value) {
  #x[[i]] <- value
  slot(object = x, name = i) = value
  return(x)
}


#' @export
#' @method [[ Assay
#'
"[[.Assay" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- 'name'
  }
  data <- slot(object = x, name = i)
  return(data)
}

#' @method show- Assay
setMethod(
  f = 'show',
  signature = 'Assay',
  definition = function(object) {
    cat(sprintf('Assay: %s\n',object@assay_name))
    cat(sprintf('Data:\n'))
    cat(str(object@data_layers,max.level = 2))
    cat(sprintf('Analysis:\n'))
    cat(str(object@analysis_layers,max.level = 2))
    
  }
)



#' Create Assay Object
#'
#' @param cell_annotations 
#' @param feature_annotations 
#' @param assay_name assay name (dna, cnv, protein)
#' @param cell_annotations table of data annotating cells, including sample labels
#' @param feature_annotations table of data annotating features. i.e. variant names or amplicon name
#'
#' @return Assay object
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' genotypes <- extract_genotypes(loom, barcodes, gt.mask=TRUE)
#' assay <- create_assay(genotypes, "snv", "sample")
#' }
create_assay<- function(assay_name, cell_annotations, feature_annotations) {
  
  
  assay <- methods::new(Class = 'Assay',
                        assay_name = assay_name,
                        cell_annotations = cell_annotations,
                        feature_annotations = feature_annotations
  )
  return(assay)
}


#' Add additional layers of data to Assay Object
#'
#' @param assay Assay object to add data to
#' @param layer_name name of layer
#' @param data new data to add
#'
#' @return Assay Object
#' @export
#'
add_data_layer<- function(assay, layer_name, data) {
  
  current_feature_names = as.character(assay@feature_annotations$id)
  new_feature_names = as.character(colnames(data))
  
  if (!all.equal(current_feature_names, new_feature_names)) {
    stop('New feature does not have same column names.')
  }
  
  if(!all.equal(dim(data), dim(assay))){
    stop(paste0('Annotations not the same length as features.'))
  }
  suppressWarnings(rownames(data) <- assay@cell_annotations$barcode)
  assay@data_layers[[layer_name]] = as_tibble(data, rownames=NA)
  
  return(assay)
}


#' Add additional layers of analysis to Assay Object. The function will check the new data has same number of cells
#'
#' @param assay Assay object to add data to
#' @param layer_name name of layer
#' @param data new data to add
#'
#' @return Assay Object
#' @export
#'
add_analysis_layer<- function(assay, layer_name, data) {
  
  data = as_tibble(data)
  
  if(nrow(assay) != nrow(data)){
    stop(paste0("analysis layer must have same number or rows (cells) as assay."))
  }
  
  suppressWarnings(rownames(data) <- assay@cell_annotations$barcode)
  assay@analysis_layers[[layer_name]] = as_tibble(data, rownames=NA)
  
  return(assay)
}

