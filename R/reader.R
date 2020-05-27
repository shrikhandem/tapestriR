
ASSAY_NAME_VARIANT = 'dna'


#' Read loom file created by Tapestri pipeline
#'
#' @param filename loom file path
#' @param min_mutation_rate only read variants with mutations rate higher then treshold. Primarly done to reduce size of data in memory
#' @return Tapestri object
#' @export
#' @examples
#' \dontrun{
#' tapestri_raw = loom_reader(filename,min_mutation_rate = 0.1)
#' }
read_loom <- function(filename, min_mutation_rate=0.05) {
  
  # filename <- "data/PE11.cells.loom"
  # min_mutation_rate = 0.1
  # experiment_name = basename(filename)
  
  h5f = rhdf5::H5Fopen(filename)
  
  features=as_tibble(h5f$row_attrs)
  cells = h5f$col_attrs
  
  
  layers = h5f$layers
  ngt = h5f$matrix
  layers$NGT = ngt
  cell_annotations = tibble(barcode=1:nrow(ngt))
  
  
  mutation_rate = apply(ngt,2, FUN = function(x) {
    sum(x %in% c(1,2)) / ncol(ngt)
  })
  mutated_variants = mutation_rate > min_mutation_rate
  filtered_features = features[mutated_variants,]
  
  assay_name = ASSAY_NAME_VARIANT
  
  assay = create_assay(assay_name = assay_name,
                       cell_annotations = cell_annotations,
                       feature_annotations = filtered_features)
  
  for(layer in names(layers)) {
    data = layers[[layer]][,mutated_variants]
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay,layer_name = layer, data = data)
  }
  
  
  rhdf5::H5Fclose(h5f)
  return(assay)
}


#' Read Tapestri multi-omics h5 file
#'
#' @param filename h5 file path
#' @param assay_name name of assay to load from h5
#' @param min_mutation_rate only read variants with mutations rate higher then treshold. Primarly done to reduce size of data in memory
#'
#' @return Tapestri multi-omics object
#' @export
#' @examples
#' \dontrun{
#' tapestri_raw = h5_reader(filename,min_mutation_rate = 0.1)
#' }
read_h5 <- function(filename, assay_name, min_mutation_rate = 0.01) {
  
  # filename <- "~/Google Drive/launches/r_package/insights_v3/data/ABseq021.h5"
  # assay_name=ASSAY_NAME_VARIANT
  # layer = 'NGT'
  # min_mutation_rate = 0.1

  h5f = rhdf5::H5Fopen(filename)
  #cell_annotations <- h5f$assays$dna$ra
  
  if(assay_name == ASSAY_NAME_VARIANT) {
    
    # filter data to make it managable in R
    ngt = h5f$assays$dna$layers$NGT
    mutation_rate = apply(ngt,1, FUN = function(x) {
      sum(x %in% c(1,2)) / ncol(ngt)
    })
    filters = mutation_rate > min_mutation_rate
  } else {
    filters = TRUE
  }
  
  h5_assay = h5f$assays[[assay_name]]
  
  filtered_features <- as_tibble(h5_assay$ca)[filters,]
  cell_annotations <- as_tibble(h5_assay$ra) %>% mutate_all(as.character)
  
  assay = create_assay(assay_name = assay_name,
                     cell_annotations = cell_annotations,
                     feature_annotations = filtered_features)

  
  for(layer in names(h5_assay$layers)) {
    data = t(h5_assay$layers[[layer]][filters, ])
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay, layer_name = layer, data = data)
  }

  rhdf5::H5Fclose(h5f)
  devnull <- base::gc()
  return(assay)     
}
