

#' The Tapestri_Multiomics Class
#'
#' Store Tapestri_Multiomics-multisample data generated from various Tapestri pipelines
#'
#' @slot cell_annotations cell annotations including sample labels
#' @slot assays 
#' @slot metadata 
#'
#' @name Tapestri_Multiomics-class
#' @rdname Tapestri_Multiomics-class
#' @exportClass Tapestri_Multiomics
#'
Tapestri_Multiomics <- setClass(
  Class = "Tapestri_Multiomics",
  slots = c(
    metadata ="list",
    assays = "list",
    cell_annotations = "tbl_df" 
  )
  #,contains = c()
)



#' @export
#' @import methods
#' @method dim Tapestri_Multiomics
#'
dim.Tapestri_Multiomics <- function(x) {
  c(nrow(x@cell_annotations))
}

#' get names of object
#' 
#' @param x 
#'
#' @export
#' @import methods
#' @method names Tapestri_Multiomics
#'
names.Tapestri_Multiomics <- function(x) {
  slotNames(x)
}

#
#' @export
#'
#' @import utils 
#' @import methods
.DollarNames.Tapestri_Multiomics <- function(x, pattern = "") {
  grep(pattern, slotNames(x), value=TRUE)
}

#
#' @export
#' @import methods
#' 
"$.Tapestri_Multiomics" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#' @export
#'
#' @import methods
"$<-.Tapestri_Multiomics" <- function(x, i, ..., value) {
  slot(object = x, name = i) = value
  return(x)
}


#' @export
#' @method [[ Tapestri_Multiomics
#' @import methods
#' 
"[[.Tapestri_Multiomics" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- 'name'
  }
  data <- slot(object = x, name = i)
  return(data)
}

#' @export
#' @method [[<- Tapestri_Multiomics
#' @import methods
#' 
"[[<-.Tapestri_Multiomics" <- function(x, i, ..., value) {
  slot(object = x, name = i) = value
  return(x)
}

#' @method show- Tapestri_Multiomics
setMethod(
  f = 'show',
  signature = 'Tapestri_Multiomics',
  definition = function(object) {
    cat(sprintf('Experiment: %s\n', object@metadata$experiment_name))
    cat(sprintf('Num Cells: %s\n', nrow(object@cell_annotations)))
    for(assay_name in names(object@assays)) {
      str(object@assays[[assay_name]], max.level = 2,give.attr = FALSE)
    }
  }
)


#' Create Multiomics Object
#'
#' @param experiment_name Name of Experiment
#' @param cell_annotations table of cell annotations. All in this object must match these cell annotations
#'
#' @return Tapestri_Multiomics Object (moo)
#' @importFrom methods new
#' @export
create_moo<- function(experiment_name, cell_annotations) {
  
  metadata = list(experiment_name = experiment_name)
  metadata[['cell_info']] = cell_annotations %>% group_by(sample) %>% summarise(cells=n())
  moo <- methods::new(Class = 'Tapestri_Multiomics',
                        metadata = metadata,
                        cell_annotations = cell_annotations
                      )
  return(moo)
}


#' Add Assay Object to Multiomics object (moo)
#'
#' @param moo Multiomics object 
#' @param assay new Assay object to add
#' @param keep_common_cells default FALSE. Will throw error if number of cells and cell barcodes dont match. If TRUE merge assays and only keep intersect of cells from all assays. 
#'
#' @return
#' @export
#'
add_assay <- function(moo, assay, keep_common_cells=FALSE) {
  
  #moo = experiment 
  #assay = cnv
  
  if(class(assay) != "Tapestri_Assay") {
    stop('Not a valid assay.')  
  }

  #check if the cells in the new assay matches current assay cells
  if(length(moo@cell_annotations$barcode) != length(assay@cell_annotations$barcode) ||
    all.equal(moo@cell_annotations$barcode, assay@cell_annotations$barcode) != TRUE) {
    
    #if not, but we want to keep only the intersect
    if (keep_common_cells) {
      interected_barcodes = dplyr::intersect(assay@cell_annotations$barcode, moo@cell_annotations$barcode)
      
      suppressWarnings(rownames(moo@cell_annotations) <- moo@cell_annotations$barcode)
      moo@cell_annotations = moo@cell_annotations[interected_barcodes,]
      
      #new_moo = create_moo(experiment_name = moo@experiment_name,cell_annotations = interected_cell_annotations)

      for(assay_name in names(moo@assays)) {
        a = cell_filter_assay(assay = moo@assays[[assay_name]], barcodes = interected_barcodes)
        moo@assays[[assay_name]] = a
      }
      
      moo@assays[[assay@assay_name]] = cell_filter_assay(assay = assay, barcodes = interected_barcodes)
    
    } else {
      stop('Cell IDs do not match')  
    }
  } else {
    #the cell info is the same between current assays and new one
    moo@assays[[assay@assay_name]] = assay    
    
  }
  moo@metadata[['cell_info']] = moo@cell_annotations %>% group_by(sample) %>% summarise(cells=n())
  

  return(moo) 
}


cell_filter_assay<- function(assay, barcodes) {
  
  for(layer in names(assay@data_layers)){
    suppressWarnings(rownames(assay@data_layers[[layer]]) <- assay@cell_annotations$barcode)
    assay@data_layers[[layer]] = assay@data_layers[[layer]][barcodes,]
  }
  for(layer in names(assay@analysis_layers)){
    suppressWarnings(rownames(assay@analysis_layers[[layer]]) <- assay@cell_annotations$barcode)
    assay@analysis_layers[[layer]] = assay@analysis_layers[[layer]][barcodes,]
  }
  suppressWarnings(rownames(assay@cell_annotations) <- assay@cell_annotations$barcode)
  assay@cell_annotations = assay@cell_annotations[barcodes,]
  
  assay@metadata[['cell_info']] = assay@cell_annotations %>% group_by(sample) %>% summarise(cell=n())
  
  return(assay)
}



