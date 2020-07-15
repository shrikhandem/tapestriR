
ASSAY_NAME_VARIANT = 'dna'


#' Read Loom file created by Tapestri Pipeline
#'
#' @param filename Loom file path
#' @param min_mutation_rate only read variants with mutation rates higher than the threshold, primarily done to reduce the size of data in memory
#' @return Tapestri Object
#' @export
#' @import rhdf5
#' 
read_loom <- function(filename, min_mutation_rate=0.05) {
  
  # filename <- "data/PE11.cells.loom"
  # min_mutation_rate = 0.1
  # experiment_name = basename(filename)
  
  assay_name = ASSAY_NAME_VARIANT
  h5f = rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY")
  

  object_str = rhdf5::h5ls(h5f)
  object_str$object_names = paste0(object_str$group,'/',object_str$name)
  dims = stringr::str_match(string = object_str$dim, pattern = '(.*) x (.*)')
  object_str$rows = as.numeric(dims[,2])
  object_str$columns = as.numeric(dims[,3])
  ngt_path = 'matrix'
  nvariants = object_str$columns[object_str$object_names == ngt_path | object_str$name == ngt_path][1]  
  mutation_rate = c()
  chunks = c(seq(from=1, to=nvariants,by = 10000), nvariants)
  for(i in 2:length(chunks)){
    ngt = rhdf5::h5read(h5f,ngt_path, index = list(NULL,chunks[i-1]:chunks[i]))      
    mutated <- (ngt == 1 | ngt == 2)
    mutation_rate <- c(mutation_rate,
                       base::colMeans(mutated, na.rm = T))
    if(!all.equal(c(0,1,2,3), sort(unique(c(ngt))))) stop('NGT layer can only contain 0, 1, 2, 3 values.')  
  }
  
  mutated_variants = which(mutation_rate > min_mutation_rate)
  

  cell_annotations <- as_tibble(rhdf5::h5read(h5f,sprintf("col_attrs")))
  
  if(nrow(cell_annotations) ==0) {
    cell_annotations = tibble(
      sample=rep(basename(filename),nrow(ngt)),
      barcode=as.character(1:nrow(ngt))
    )
  }
  cell_annotations <- validate_cell_annotations(cell_annotations = cell_annotations)  

  
  all_features =as_tibble(rhdf5::h5read(h5f,sprintf("row_attrs")))
  filtered_features = all_features[mutated_variants,]
  filtered_features = validate_variant_annotations(filtered_features)
  
  assay = create_assay(assay_name = assay_name,
                       cell_annotations = cell_annotations,
                       feature_annotations = filtered_features)
  
  filtered_ngt = data = rhdf5::h5read(h5f,ngt_path, index = list(NULL,mutated_variants))
  colnames(filtered_ngt) <- filtered_features$id
  assay = add_data_layer(assay = assay,layer_name = 'NGT', data = filtered_ngt)
  
  layer_names = rhdf5::h5ls(h5f&sprintf("layers"))
  for(layer in layer_names$name) {
    
    data = rhdf5::h5read(h5f,sprintf("layers/%s", layer), index = list(NULL,mutated_variants))
    
    #data = layers[[layer]][,mutated_variants]
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay, layer_name = layer, data = data)
  }
  
  
  rhdf5::h5closeAll()
  return(assay)
}

#' Read Tapestri Multi-omics H5 file
#'
#' @param filename H5 file path
#' @param assay_name name of assay to load from H5
#' @param min_mutation_rate only read variants with mutation rates higher than the threshold, primarily done to reduce the size of data in memory
#'
#' @return Tapestri Multiomics Object
#' @export
#' @import rhdf5
#' @import stringr
#' @examples
#' \dontrun{
#' tapestri_raw = h5_reader(filename,min_mutation_rate = 0.1)
#' }
read_assay_h5 <- function(filename, assay_name, min_mutation_rate = 0.005) {
  
  #filename <- "~/Google Drive/launches/r_package/data/merged_all.h5"
  #assay_name='dna'
  # layer = 'NGT'
  # min_mutation_rate = 0.1

  h5f = rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY")
  object_str = rhdf5::h5ls(h5f)
  object_str$object_names = paste0(object_str$group,'/',object_str$name)
  dims = stringr::str_match(string = object_str$dim, pattern = '(.*) x (.*)')
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
  
  cell_annotations <- as_tibble(rhdf5::h5read(h5f,ra))
  cell_annotations <- validate_cell_annotations(cell_annotations = cell_annotations)  
  
  features =as_tibble(rhdf5::h5read(h5f,ca))
  
  if(assay_name == ASSAY_NAME_VARIANT) {
    # filter data to make it manageable in R
    
    features = validate_variant_annotations(features)  
    
    ngt_path = sprintf("/assays/%s/layers/NGT",assay_name)    
    if(!ngt_path %in% object_str$object_names) {
      stop(sprintf('H5 file does not seem to be a standard structure. Missing %s.', ngt_path))
    }
    
    nvariants = object_str$rows[object_str$object_names == ngt_path]
    mutation_rate = c()
    chunks = c(seq(from=1, to=nvariants,by = 10000),nvariants)
    
    for(i in 2:length(chunks)){
      ngt = rhdf5::h5read(h5f,ngt_path, index = list(chunks[i-1]:chunks[i],NULL))
      mutated <- (ngt == 1 | ngt == 2)
      mutation_rate <- c(mutation_rate,
                         base::rowMeans(mutated, na.rm = T))
      if(!all.equal(c(0,1,2,3), sort(unique(c(ngt))))) stop('NGT layer can only contain 0, 1, 2, 3 values.')
    }
    filters = which(mutation_rate > min_mutation_rate)
    
    
  } else {
    
    filters = 1:nrow(features)
  }
  

  filtered_features <- features[filters,]
  
  assay = create_assay(assay_name = assay_name,
                     cell_annotations = cell_annotations,
                     feature_annotations = filtered_features)

  layer_names = rhdf5::h5ls(h5f&sprintf("assays/%s/layers",assay_name))
  
  for(layer in layer_names$name) {
    
    data = rhdf5::h5read(h5f,sprintf("assays/%s/layers/%s",assay_name,layer), index = list(filters,NULL))
    
    data = t(data)
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay, layer_name = layer, data = data)
  }
  
  rhdf5::h5closeAll()
  
  return(assay)     
}


#' Read a TAP file generated by Tapestri Pipeline
#'
#' @param filename path to TAP file
#' @param experiment_name name of experiment, if none file name is used
#'
#' @return Tapestri Multiomics Object
#' @export
#'
read_tap <- function(filename, experiment_name = NA) {
  
  if(is.na(experiment_name)) experiment_name = basename(filename)
  
  h5f = rhdf5::H5Fopen(filename, flags = "H5F_ACC_RDONLY")
  object_str = rhdf5::h5ls(h5f)
  object_str$object_names = paste0(object_str$group,'/',object_str$name)
  dims = stringr::str_match(string = object_str$dim, pattern = '(.*) x (.*)')
  object_str$rows = as.numeric(dims[,2])
  object_str$columns = as.numeric(dims[,3])
  assay_names = rhdf5::h5ls(h5f&sprintf("/assays"), recursive = FALSE)
  
  a = read_assay_h5(filename = filename, assay_name = assay_names$name[1])
  moo = create_moo(experiment_name = experiment_name, cell_annotations = a@cell_annotations)
  moo = add_assay(moo,a) 
  for(assay in assay_names$name[-1]) {
    a = read_assay_h5(filename = filename, assay_name = assay)
    moo = add_assay(moo,a) 
  }
  return(moo)
}



#' Read variant assay from data exported from Tapestri Insights
#'
#' Creates an Assay Object from the Tapestri Insights exported zip file. Includes all the DNA variant layers plus variant annotations. If subclone information exists, creates an analysis layer with the subclone labels.
#'
#' @param export_dir directory path to exported data. Should contain "AF.csv", "DP.csv", "GQ.csv", "NGT.csv", "README.txt", "Variants.csv"
#'
#' @return
#' @export
#'
read_insights_export <- function(export_dir) {
  layers = c('NGT', 'AF', 'DP', 'GQ')
  annotations = 'Variants'
  files = c(layers, annotations)
  files_exist = file.exists(paste0(export_dir, '/', files, '.csv'))
  if (any(!files_exist))
    stop(sprintf('Some files missing: %s', paste(files[!files_exist])))
  
  
  ngt = readr::read_csv(paste0(export_dir, '/', 'NGT.csv'))

  # check if NGT file has only relavent values
  #unique(c(as.matrix(ngt %>% select(-Sample, -Cell))))
  if(!all.equal(c(0,1,2,3), sort(unique(c(as.matrix(ngt %>% select(-matches(c('Sample','Cell','Subclone'))))))))) stop('NGT layer can only contain 0, 1, 2, 3 values.')
  
  cell_annotations = ngt %>% select(sample = Sample, barcode = Cell) %>%
    mutate(barcode = as.character(barcode),
           id = paste0(sample, '-', barcode))
  
  cell_annotations <- validate_cell_annotations(cell_annotations = cell_annotations)  
  
  variant_annotations = readr::read_csv(paste0(export_dir, '/', annotations, '.csv'))
  variant_annotations = variant_annotations %>% mutate(id=variant_annotations$Variant)
  variant_annotations = validate_variant_annotations(variant_annotations)
  
  
  variants_from_insights = create_assay(
    assay_name = 'variants',
    cell_annotations = cell_annotations,
    feature_annotations = variant_annotations
  )
  suppressMessages(for (layer in layers) {
    layer_data = readr::read_csv(paste0(export_dir, '/', layer, '.csv'))
    layer_data = layer_data %>% select(-Sample,-Cell) %>% select(variant_annotations$id)
    
    variants_from_insights = add_data_layer(assay = variants_from_insights,
                                            layer_name = layer,
                                            data = layer_data)
    
  })
  
  if ('Subclone' %in% colnames(ngt)) {
    variants_from_insights = add_analysis_layer(assay=variants_from_insights,layer_name = 'Subclone',data = ngt$Subclone)
  }
  
  return(variants_from_insights)  
}



validate_cell_annotations <- function(cell_annotations){
  
  cell_annotations = as_tibble(cell_annotations) %>% mutate_all(as.character)
  
  if(!'sample' %in% colnames(cell_annotations)){
    cell_annotations = cell_annotations %>% mutate(sample=basename(filename))
  }
  if(!'barcode' %in% colnames(cell_annotations)){
    cell_annotations = cell_annotations %>% mutate(barcode=row_number())
  }

  if(!'id' %in% colnames(cell_annotations)){
    cell_annotations = cell_annotations %>% mutate(id=paste0(sample,'-',barcode))
  } else if (any(duplicated(cell_annotations$id))) {
    stop('Cell annotations id must be unique.')
  }

  return(cell_annotations)
}


#' Annotate variants based on variant id field
#'
#' @param variant_annotations table of annotations, unique id column must exist, if it is in format 'SF3B1:chr2:198266943:C/T', try to parse into CHROM and POS
#'
#' @return
#' @export
#' @import stringr
#' @import tibble

validate_variant_annotations <- function(variant_annotations) {
  
  if(!'id' %in% colnames(variant_annotations)) stop('id column must exist.')
  if(any(duplicated(variant_annotations$id))) {
    warning("Duplicate variant ids present. You're likely running an older experiment. 
            Don't worry. We'll fix the duplicate issue here, and new runs will not have this issue.")
    
    variant_annotations = variant_annotations %>% distinct(id, .keep_all = TRUE)
    }
  
  if(any(!c('CHROM','POS') %in% colnames(variant_annotations))){
    split_ids = stringr::str_match(variant_annotations$id,'(.*?):?chr(.*?):([[:digit:]]*):')
    if(any(is.na(split_ids[,3])) | any(is.na(split_ids[,4]))) stop(sprintf('id column cannot be easily parsed into CHROM and POS.'))
       
    variant_annotations = variant_annotations %>% mutate(CHROM=split_ids[,3], POS=as.numeric(split_ids[,4]), gene_name=split_ids[,2])  
    variant_annotations = variant_annotations %>% arrange(as.numeric(CHROM), as.numeric(POS))
  }

  return(variant_annotations)
}

