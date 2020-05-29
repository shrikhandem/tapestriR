
ASSAY_NAME_VARIANT = 'dna'


#' Read loom file created by Tapestri pipeline
#'
#' @param filename loom file path
#' @param min_mutation_rate only read variants with mutations rate higher then threshold. Primary done to reduce size of data in memory
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
  
  assay_name = ASSAY_NAME_VARIANT
  h5f = rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY")
  #h5ls(h5f)
  
  features =as_tibble(h5read(h5f,sprintf("row_attrs")))
  ngt = h5f$matrix
  
  cell_annotations <- as_tibble(h5read(h5f,sprintf("col_attrs"))) %>% mutate_all(as.character)
  if(nrow(cell_annotations) ==0) {
    cell_annotations = tibble(
      sample=rep(basename(filename),nrow(ngt)),
      barcode=as.character(1:nrow(ngt))
    )
  }

  # mutation_rate = apply(ngt,2, FUN = function(x) {
  #   sum(x %in% c(1,2)) / ncol(ngt)
  # })
   
  mutated <- (ngt == 1 | ngt == 2)
  mutation_rate <- base::colMeans(mutated, na.rm = T)

  mutated_variants = which(mutation_rate > min_mutation_rate)
  
  
  filtered_features = features[mutated_variants,]
  
  assay = create_assay(assay_name = assay_name,
                       cell_annotations = cell_annotations,
                       feature_annotations = filtered_features)
  
  filtered_ngt = ngt[,mutated_variants]
  colnames(filtered_ngt) <- filtered_features$id
  assay = add_data_layer(assay = assay,layer_name = 'NGT', data = filtered_ngt)
  
  layer_names = h5ls(h5f&sprintf("layers"))
  for(layer in layer_names$name) {
    
    data = h5read(h5f,sprintf("layers/%s", layer), index = list(NULL,mutated_variants))
    
    #data = layers[[layer]][,mutated_variants]
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay,layer_name = layer, data = data)
  }
  
  
  h5closeAll()
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
#' @import rhdf5
#' @examples
#' \dontrun{
#' tapestri_raw = h5_reader(filename,min_mutation_rate = 0.1)
#' }
read_h5 <- function(filename, assay_name, min_mutation_rate = 0.01) {
  
  #filename <- "~/Google Drive/launches/r_package/data/merged_all.h5"
  # assay_name=ASSAY_NAME_VARIANT
  # layer = 'NGT'
  # min_mutation_rate = 0.1

  h5f = rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY")
  #h5ls(h5f)
  
  features =as_tibble(h5read(h5f,sprintf("assays/%s/ca",assay_name)))
  cell_annotations <- as_tibble(h5read(h5f,sprintf("assays/%s/ra",assay_name))) %>% mutate_all(as.character)
  if(!'sample' %in% colnames(cell_annotations)){
    cell_annotations = cell_annotations %>% mutate(sample=basename(filename))
  }
  
  if(assay_name == ASSAY_NAME_VARIANT) {
    
    # filter data to make it manageable in R
    
    #ngt = h5f$assays$dna$layers$NGT
    #ngt = h5f$"assays/dna/layers/NGT"
    
    ngt = h5read(h5f,"assays/dna/layers/NGT")

    #### todo: faster way to filter large data
    # mutation_rate = apply(ngt,1, FUN = function(x) {
    #   sum(x %in% c(1,2)) / ncol(ngt)
    # })
    # 
    #mutated <- (ngt == 1 | ngt == 2)
    #rowMeans(mutated) * 100
    
    mutated <- (ngt == 1 | ngt == 2)
    mutation_rate <- base::rowMeans(mutated, na.rm = T)
    filters = which(mutation_rate > min_mutation_rate)
    
  } else {
    
    filters = 1:nrow(features)
  }
  

  filtered_features <- features[filters,]
  
  assay = create_assay(assay_name = assay_name,
                     cell_annotations = cell_annotations,
                     feature_annotations = filtered_features)

  layer_names = h5ls(h5f&sprintf("assays/%s/layers",assay_name))
  
  for(layer in layer_names$name) {
    
    data = h5read(h5f,sprintf("assays/%s/layers/%s",assay_name,layer), index = list(filters,NULL))
    
    data = t(data)
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay, layer_name = layer, data = data)
  }

  
  h5closeAll()
  
  return(assay)     
}
