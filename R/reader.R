
ASSAY_NAME_VARIANT = 'dna'


#' Read loom file created by Tapestri pipeline
#'
#' @param filename loom file path
#' @param min_mutation_rate only read variants with mutations rate higher then threshold. Primary done to reduce size of data in memory
#' @return Tapestri object
#' @export
#' @import rhdf5
#' @examples
#' \dontrun{
#' tapestri_raw = loom_reader(filename, min_mutation_rate = 0.1)
#' }
read_loom <- function(filename, min_mutation_rate=0.05) {
  
  # filename <- "data/PE11.cells.loom"
  # min_mutation_rate = 0.1
  # experiment_name = basename(filename)
  
  assay_name = ASSAY_NAME_VARIANT
  h5f = rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY")
  #rhdf5::h5ls(h5f)
  
  ngt = h5f$matrix
  mutated <- (ngt == 1 | ngt == 2)
  mutation_rate <- base::colMeans(mutated, na.rm = T)
  mutated_variants = which(mutation_rate > min_mutation_rate)
  # mutation_rate = apply(ngt,2, FUN = function(x) {
  #   sum(x %in% c(1,2)) / ncol(ngt)
  # })
  
  features =as_tibble(rhdf5::h5read(h5f,sprintf("row_attrs")))
    
  cell_annotations <- as_tibble(rhdf5::h5read(h5f,sprintf("col_attrs")))

  if(nrow(cell_annotations) ==0) {
    cell_annotations = tibble(
      sample=rep(basename(filename),nrow(ngt)),
      barcode=as.character(1:nrow(ngt))
    )
  }
  cell_annotations <- validate_cell_annotations(cell_annotations = cell_annotations)  


  filtered_features = features[mutated_variants,]
  
  assay = create_assay(assay_name = assay_name,
                       cell_annotations = cell_annotations,
                       feature_annotations = filtered_features)
  
  filtered_ngt = ngt[,mutated_variants]
  colnames(filtered_ngt) <- filtered_features$id
  assay = add_data_layer(assay = assay,layer_name = 'NGT', data = filtered_ngt)
  
  layer_names = rhdf5::h5ls(h5f&sprintf("layers"))
  for(layer in layer_names$name) {
    
    data = rhdf5::h5read(h5f,sprintf("layers/%s", layer), index = list(NULL,mutated_variants))
    
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
#' @import tidyverse
#' @examples
#' \dontrun{
#' tapestri_raw = h5_reader(filename,min_mutation_rate = 0.1)
#' }
read_h5 <- function(filename, assay_name, min_mutation_rate = 0.01) {
  
  #filename <- "~/Google Drive/launches/r_package/data/merged_all.h5"
  #assay_name=ASSAY_NAME_PROTEIN
  # layer = 'NGT'
  # min_mutation_rate = 0.1

  h5f = rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY")
  object_str = h5ls(h5f)
  object_str$object_names = paste0(object_str$group,'/',object_str$name)
  dims = str_match(string = object_str$dim, pattern = '(.*) x (.*)')
  object_str$rows = as.numeric(dims[,2])
  object_str$columns = as.numeric(dims[,3])
  
  ca = sprintf("/assays/%s/ca",assay_name)
  ra = sprintf("/assays/%s/ra",assay_name)
  
  
  if(!ca %in% object_str$object_names) {
    stop(sprintf('H5 file does not seem to be a standard structure. Missing %s.', ca))
  }
  if(!ra %in% object_str$object_names) {
    stop(sprintf('H5 file does not seem to be a standard structure. Missing %s.', ra))
  }
  
  features =as_tibble(rhdf5::h5read(h5f,ca))
  cell_annotations <- as_tibble(rhdf5::h5read(h5f,ra))
  
  cell_annotations <- validate_cell_annotations(cell_annotations = cell_annotations)  
  
  if(assay_name == ASSAY_NAME_VARIANT) {
    
    # filter data to make it manageable in R
    ngt_path = sprintf("/assays/%s/layers/NGT",assay_name)    
    if(!ngt_path %in% object_str$object_names) {
      stop(sprintf('H5 file does not seem to be a standard structure. Missing %s.', ngt_path))
    }
    
    mutation_rate = c()
    chunks = c(seq(from=1, to=object_str$rows[object_str$object_names == ngt_path],by = 10000), object_str$rows[object_str$object_names == ngt_path])
    for(i in 2:length(chunks)){
      ngt = rhdf5::h5read(h5f,ngt_path, index = list(chunks[i-1]:chunks[i],NULL))      
      mutated <- (ngt == 1 | ngt == 2)
      mutation_rate <- c(mutation_rate,
                         base::rowMeans(mutated, na.rm = T))
      
    }
    filters = which(mutation_rate > min_mutation_rate)

    #### todo: faster way to filter large data
    # mutation_rate = apply(ngt,1, FUN = function(x) {
    #   sum(x %in% c(1,2)) / ncol(ngt)
    # })
    # 
    #mutated <- (ngt == 1 | ngt == 2)
    #rowMeans(mutated) * 100
    
    
  } else {
    
    filters = 1:nrow(features)
  }
  

  filtered_features <- features[filters,]
  
  assay = create_assay(assay_name = assay_name,
                     cell_annotations = cell_annotations,
                     feature_annotations = filtered_features)

  layer_names = h5ls(h5f&sprintf("assays/%s/layers",assay_name))
  
  for(layer in layer_names$name) {
    
    data = rhdf5::h5read(h5f,sprintf("assays/%s/layers/%s",assay_name,layer), index = list(filters,NULL))
    
    data = t(data)
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay, layer_name = layer, data = data)
  }

  
  h5closeAll()
  
  return(assay)     
}

validate_cell_annotations <- function(cell_annotations){
  
  cell_annotations = as_tibble(cell_annotations) %>% mutate_all(as.character)
  
  if(!'sample' %in% colnames(cell_annotations)){
    cell_annotations = cell_annotations %>% mutate(sample=basename(filename))
  }
  if(!'barcode' %in% colnames(cell_annotations)){
    cell_annotations = cell_annotations %>% mutate(barcode=row_number())
  }
  return(cell_annotations)
}


#' Annotate Variants based on variant ID field
#'
#' @param variant_ids variant id like 'ASXL1:chr20:31022959:T/C' 
#'
#' @return
#' @export
#' @import tidyverse

annotate_variants <- function(variant_ids) {
  
  variant_annotations = str_match(variant_ids,'(.*?):?chr(.*?):([[:digit:]]*):')
  variant_annotations = tibble(id=variant_ids, CHROM=variant_annotations[,3], POS=as.numeric(variant_annotations[,4]), gene_name=variant_annotations[,2])
  
  return(variant_annotations)
}


