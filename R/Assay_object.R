
setOldClass(c("tbl_df", "tbl", "data.frame"))

#' Tapestri_Assay class
#'
#' Store assay specific details for given sample
#'
#' @slot assay_name 
#' @slot cell_annotations 
#' @slot feature_annotations 
#' @slot data_layers 
#' @slot analysis_layers 
#' @slot metadata 
#'
#' @name Tapestri_Assay-class
#' @rdname Tapestri_Assay-class
#' @import dplyr
#' @exportClass Tapestri_Assay
Tapestri_Assay <- setClass(
  Class = "Tapestri_Assay",
  slots = c(
    assay_name = "character",
    metadata = 'list',
    feature_annotations = "tbl_df",
    cell_annotations = "tbl_df",
    data_layers = "list",
    analysis_layers = "list"
  )
  #,contains = c(class(tibble()))
)
setClassUnion(name="AssayorNULL", c("Tapestri_Assay", "NULL"))

#' @export
#' @method dim Tapestri_Assay
#'
dim.Tapestri_Assay <- function(x) {
  c(nrow(x@cell_annotations), nrow(x@feature_annotations))
}

#' @export
#' @method nrow Tapestri_Assay
#'
nrow.Tapestri_Assay <- function(x) {
  nrow(x@cell_annotations)
}


#' @export
#'
"$.Tapestri_Assay" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#
#' @export
#' @import utils
.DollarNames.Tapestri_Assay <- function(x, pattern = "") {
  grep(pattern, slotNames(x), value=TRUE)
}


#' @export
#'
"$<-.Tapestri_Assay" <- function(x, i, ..., value) {
  #x[[i]] <- value
  slot(object = x, name = i) = value
  return(x)
}


#' @export
#' @method [[ Tapestri_Assay
#'
"[[.Tapestri_Assay" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- 'name'
  }
  data <- slot(object = x, name = i)
  return(data)
}


#' @export
#' @method [[<- Tapestri_Assay
#'
"[[<-.Tapestri_Assay" <- function(x, i, ..., value) {
  slot(object = x, name = i) = value
  return(x)
}


#' @method show- Tapestri_Assay
setMethod(
  f = 'show',
  signature = 'Tapestri_Assay',
  definition = function(object) {
    cat(str(object,max.level = 3, give.attr = FALSE, vec.len = 2))

  }
)



#' Create Tapestri_Assay Object
#'
#' @param cell_annotations 
#' @param feature_annotations 
#' @param assay_name assay name (dna, cnv, protein)
#' @param cell_annotations table of data annotating cells, including sample labels
#' @param feature_annotations table of data annotating features. i.e. variant names or amplicon name
#'
#' @return Tapestri_Assay object
#' @importFrom methods new
#' @export
#' 
create_assay<- function(assay_name, cell_annotations, feature_annotations) {
  
  if(!'id' %in% colnames(feature_annotations)) {
    stop('id column must exist in feature_annotations.')
  }
  
  if(!'sample' %in% colnames(cell_annotations)) {
    stop('sample column must exist cell_annotations')
  }

  if(!'barcode' %in% colnames(cell_annotations)) {
    stop('barcode column must exist cell_annotations')
  }
  metadata = list()
  metadata[['cell_info']] = cell_annotations %>% group_by(sample) %>% summarise(cells=n())
  
  assay <- methods::new(Class = 'Tapestri_Assay',
                        assay_name = assay_name,
                        metadata = metadata,
                        cell_annotations = cell_annotations,
                        feature_annotations = feature_annotations
  )
  return(assay)
}


#' Add additional layers of data to Tapestri_Assay Object
#'
#' @param assay Tapestri_Assay object to add data to
#' @param layer_name name of layer
#' @param data new data to add
#'
#' @return Tapestri_Assay Object
#' @export
#'
add_data_layer<- function(assay, layer_name, data) {
  
  current_feature_names = as.character(assay@feature_annotations$id)
  new_feature_names = as.character(colnames(data))
  col_check = all.equal(current_feature_names, new_feature_names)
  if (col_check != TRUE) {
    stop(sprintf('New feature does not have same column names.\n%s', col_check))
  }
  dim_check = all.equal(dim(data), dim(assay))
  if(dim_check !=TRUE){
    stop(sprintf('Demension of new data layer must be the same as current layers.\n%s', dim_check))
  }
  
  # suppressWarnings(
  #   rownames(data) <- paste(assay$cell_annotations$sample, assay@cell_annotations$barcode,sep = '_')
  # )
  assay@data_layers[[layer_name]] = as_tibble(data)
  
  return(assay)
}


#' Add additional layers of analysis to Tapestri_Assay Object. The function will check the new data has same number of cells
#'
#' @param assay Tapestri_Assay object to add data to
#' @param layer_name name of layer
#' @param data new data to add
#'
#' @return Tapestri_Assay Object
#' @export
#'
add_analysis_layer<- function(assay, layer_name, data) {
  
  data = as_tibble(data)
  
  if(nrow(assay) != nrow(data)){
    stop(paste0("analysis layer must have same number or rows (cells) as assay."))
  }
  
  assay@analysis_layers[[layer_name]] = data
  
  return(assay)
}


#' Subset assay
#'
#' @param assay Assay object to subset
#' @param keep_cell_ids vector of cell ids to keep
#' @param keep_feature_ids vector of feature ids to keep
#'
#' @return return subsetted assay
#' @export
#'
subset_assay<- function(assay, keep_cell_ids=TRUE, keep_feature_ids = TRUE) {
  
  if (length(keep_cell_ids) > 1) {
    cell_ind = match(keep_cell_ids, assay@cell_annotations$id)
    if (any(is.na(cell_ind))) stop('cell ids to keep dont exist in assay.')
  } else if(keep_cell_ids==TRUE) {
    cell_ind = 1:nrow(assay@cell_annotations)
  }

  if (length(keep_feature_ids) > 1) {
    feature_ind = match(keep_feature_ids, assay@feature_annotations$id)
    if (any(is.na(cell_ind))) stop('features to keep dont exist in assay.')
  } else if(keep_feature_ids==TRUE) {
    feature_ind = 1:nrow(assay@feature_annotations)
  }
  
  
  for(layer in names(assay@data_layers)){
    assay@data_layers[[layer]] = assay@data_layers[[layer]][cell_ind,feature_ind]
  }
  
  for(layer in names(assay@analysis_layers)){
    assay@analysis_layers[[layer]] = assay@analysis_layers[[layer]][cell_ind,feature_ind]
  }
  
  assay@cell_annotations = assay@cell_annotations[cell_ind,]
  assay@feature_annotations = assay@feature_annotations[feature_ind,]
  
  assay@metadata[['cell_info']] = assay@cell_annotations %>% group_by(sample) %>% summarise(cell=n())
  
  return(assay)
}
