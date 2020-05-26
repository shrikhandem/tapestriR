
#' The Multiomics Class
#'
#' Store multiomics-multisample data generated from various Tapestri pipelines
#' @slot name Analysis name
#' @slot summary Sample or Assay summary
#' @slot version Package version or Tapestri pipeline version or Loom file version
#' @slot rna RNA assay
#' @slot snv SNV assay
#' @slot dna DNA assay
#' @slot prot Protein assay
#'
#' @name Multiomics-class
#' @rdname Multiomics-class
#' @exportClass Multiomics
#'
Multiomics <- setClass(
  Class = "Multiomics",
  slots = c(
    experiment_name = "character",
    dna = "Assay",
    cnv = "Assay",
    protein = "Assay",
    cell_annotations = "tbl_df" 
  ),
  contains = c()
)



#' @export
#' @method dim Multiomics
#'
dim.Multiomics <- function(x) {
  c(nrow(x@cell_annotations), nrow(x@feature_annotations))
}


#' @export
#'
"$.Multiomics" <- function(x, i, ...) {
  return(slot(object = x, name = i))
}

#' @export
#'
"$<-.Multiomics" <- function(x, i, ..., value) {
  #x[[i]] <- value
  slot(object = x, name = i) = value
  return(x)
}



#' @export
#' @method [[ Multiomics
#'
"[[.Multiomics" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- 'name'
  }
  data <- slot(object = x, name = i)
  return(data)
}

#' @method show- Multiomics
setMethod(
  f = 'show',
  signature = 'Multiomics',
  definition = function(object) {
    cat('Multiomics Object (moo)\n')
    cat(sprintf('Experiment: %s\n', object@experiment_name))
    cat(sprintf('%s: %s \n', ASSAY_NAME_VARIANT, .hasSlot(object, ASSAY_NAME_VARIANT)))
    str(object[[ASSAY_NAME_VARIANT]]@data_layers, max.level = 2)
    cat(sprintf('%s: %s \n', ASSAY_NAME_READ_COUNT, .hasSlot(object, ASSAY_NAME_READ_COUNT)))
    str(object[[ASSAY_NAME_READ_COUNT]]@data_layers, max.level = 2)
    cat(sprintf('%s: %s \n', ASSAY_NAME_PROTEIN, .hasSlot(object, ASSAY_NAME_PROTEIN)))
    str(object[[ASSAY_NAME_PROTEIN]]@data_layers, max.level = 2)
  }
)




#' Create Multiomics Object
#'
#' @param data data
#' @param type assay type
#' @param name sample name
#' @return Analyte Analyte object
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' genotypes <- extract_genotypes(loom, barcodes, gt.mask=TRUE)
#' assay <- create_assay(genotypes, "snv", "sample")
#' }
create_moo<- function(experiment_name, cell_annotations) {
  
  moo <- methods::new(Class = 'Multiomics',
                        experiment_name = experiment_name,
                        cell_annotations = cell_annotations
                      )
  return(moo)
}


add_assay <- function(moo, assay, keep_common_cells=FALSE) {
  
  # moo = ab21
  # assay = cnv
  # 
  
  #check if the cells in the new assay matches current assay cells
  if(length(moo@cell_annotations$barcode) != length(assay@cell_annotations$barcode) ||
    all.equal(moo@cell_annotations$barcode, assay@cell_annotations$barcode) != TRUE) {
    
    #if not, but we want to keep only the intersect
    if (keep_common_cells) {
      common_barcodes = tibble(barcode = dplyr::intersect(assay@cell_annotations$barcode, moo@cell_annotations$barcode))
      
      for(layer in names(assay@data_layers)){
        assay@data_layers[[layer]] = assay@data_layers[[layer]][common_barcodes$barcode,]
      }
      for(layer in names(assay@analysis_layers)){
        assay@analysis_layers[[layer]] = assay@analysis_layers[[layer]][common_barcodes$barcode,]
      }
      
      rownames(assay@cell_annotations) <- assay@cell_annotations$barcode
      assay@cell_annotations = assay@cell_annotations[common_barcodes$barcode,]
      
      slot(moo,assay@assay_name) = assay
    } else {
      stop('Cell IDs do not match')  
    }
  } else {
    #the cell info is the same between current assays and new one
    slot(moo,assay@assay_name) = assay    
  }
  
 return(moo) 
}








