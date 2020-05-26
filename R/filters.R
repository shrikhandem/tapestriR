


#' Read loom file created by Tapestri pipeline
#'
#' @param filename loom file name
#' @param min_mutation_rate only read variants with mutations rate higher then treshold. Primarly done to reduce size of data in memory
#' @return Tapestri object
#' @export
#' @examples
#' \dontrun{
#' tapestri_raw = loom_reader(filename,min_mutation_rate = 0.1)
#' }
loom_reader <- function(filename, min_mutation_rate=0.05) {
  
  # filename <- "data/PE11.cells.loom"
  # min_mutation_rate = 0.1
  
  h5f = rhdf5::H5Fopen(filename)
  
  features=as_tibble(h5f$row_attrs)
  cells = h5f$col_attrs
 
  layers = h5f$layers
  ngt = h5f$matrix
  layers$NGT = ngt
  
  
  
  mutation_rate = apply(ngt,2, FUN = function(x) {
    sum(x %in% c(1,2)) / ncol(ngt)
  })
  mutated_variants = mutation_rate > min_mutation_rate
  filtered_features = features[mutated_variants,]
  assay_name = 'variant'
  
  tapestri_object = list()
  for(layer in names(layers)) {
    data = layers[[layer]][,mutated_variants]
    colnames(data) <- filtered_features$id
    tapestri_object = add_feature(tapestri_object = tapestri_object, assay_name = assay_name, feature_name = layer, new_data = data)
  }
  tapestri_object = add_feature_annotations( tapestri_object = tapestri_object, assay_name = assay_name, annotations = filtered_features)
    
  rhdf5::H5Fclose(h5f)
  return(tapestri_object)
}

#' Add annotations to a Tapestri object. 
#'
#' @param tapestri_object 
#' @param assay_name feature belongs to assay variants, cnv, proteins 
#' @param annotations data.frame of annotations. should be same length as features in the assay. 
#'
#' @return tapestri object
#' @export
#'
#' @examples
add_feature_annotations <- function(tapestri_object, assay_name, annotations){

  if(ncol(tapestri_object[[assay_name]][['features']][[1]]) != nrow(annotations)){
    stop(paste0('Annotations not the same length as features.'))
  }
  
  if(!is.null(tapestri_object[[assay_name]][['annotations']])){
    to_add = tibble(
      tapestri_object[[assay_name]][['annotations']],
      annotations,
    .name_repair ='unique')
  } else {
    to_add = annotations
  }
  
  tapestri_object[[assay_name]][['annotations']] = annotations
  
  return(tapestri_object)
}


#' Add features to tapestri object
#'
#' @param tapestri_object tapestri object to update
#' @param assay_name feature belongs to assay variants, cnv, proteins 
#' @param feature_name name of feature
#' @param new_data new date to add to taesptri object
#'
#' @return tapestri object
#' @export
#'
add_feature <- function(tapestri_object, assay_name, feature_name, new_data) {
  
  # tapestri_object = analysis_df
  # new_data = protein_counts_norm
  # assay_name = 'protein'
  # feature_name='normalized'
  # 
  
  # if (!assay_name %in% names(tapestri_object)) {
  #   stop(paste0("Assay ", assay_name, ' does not exist.'))
  # }

  if (!assay_name %in% names(tapestri_object)) {
    tapestri_object[[assay_name]][['features']][[feature_name]] = as_tibble(new_data) 
  } else {
    current_feature_names = as.character(colnames(tapestri_object[[assay_name]][['features']][[1]]))
    new_feature_names = as.character(colnames(new_data))
    
    if (!all.equal(current_feature_names, new_feature_names)) {
      stop('New feature does not have same column names.')
    }
    tapestri_object[[assay_name]][['features']][[feature_name]] = as_tibble(new_data) 
  }
  
  return(tapestri_object)
}


#' Add analysis to tapestri object. the analysis should have same # of cells as tapestri object
#'
#' @param tapestri_object tapestri object to update
#' @param assay_name feature belongs to assay variants, cnv, proteins 
#' @param analysis_name 
#' @param new_data 
#'
#' @return
#' @export
#'
add_analysis <- function(tapestri_object, assay_name, analysis_name, new_data) {
  
  # tapestri_object = analysis_df
  # assay_name = 'variants' 
  # analysis_name =   'clusters'
  # new_data = clusters[,-1]
  # 
  if (is.null(tapestri_object[[assay_name]]$analysis)) {
    
    tmp_data = tibble(cell_id = tapestri_object$cell_id)
    tmp_data[[analysis_name]] = tibble(new_data)  
    
    tapestri_object[[assay_name]] = tibble(
      tapestri_object[[assay_name]],
      analysis = tmp_data[,-1]
      )
                                           
  } else {
    
    tmp_data = tibble(tapestri_object[[assay_name]]$analysis)
    tmp_data[[analysis_name]] = tibble(new_data)  
    
    tapestri_object[[assay_name]]$analysis = tmp_data
    
  }
  
  return(tapestri_object)
}




#' Filter raw genotypes based on quality
#'
#' @param loom Loom file
#' @param gqc Genotype quality cutoff (default 30)
#' @param dpc Read depth cutoff (default 10)
#' @param afc Allele frequency cutoff (default 20)
#' @param mv Remove variants with < mv of known values (default 50)
#' @param mc Remove variants with < mc of known values (default 50)
#' @param mm Remove variants mutated in < mm of cells (default 1)
#' @param gt.mask mask low quality GT as missing (if GQ/DP/AF lower than cutoff, default FALSE)
#' @return Filtered genotypes
#' @export
filter_variants <- function(tapestri_object, gqc = 30, dpc = 10, afc = 20, mv = 50, mc = 50, mm = 1, gt.mask = FALSE) {
  
  # tapestri_object = tapestri_raw
  # gqc = 30
  # dpc = 10
  # afc = 20
  # mv = 50
  # mc = 50
  # mm = 1
  # gt.mask = FALSE
  # 
  
  
  data = tapestri_object$variant$features
  
  gt <- as.matrix(data$NGT)
  mask <- (gt < 3)
  mutated <- (gt == 1 | gt == 2)
  
  gq <- (data$GQ >= gqc)
  
  
  dp <- data$DP
  ad <- data$AD
  af <- matrix(100, nrow = nrow(dp), ncol = ncol(dp))
  
  af[mutated] <- ad[mutated] * 100 / dp[mutated]
  af[is.na(af)] <- 0
  af <- (af >= afc)
  dp <- (dp >= dpc)
  ngt_filter <- gq & dp & af & mask
  mv.c <- base::colMeans(ngt_filter, na.rm = T) * 100
  kept_variants <- base::which(mv.c >= mv)
  mc.c <- base::rowMeans(ngt_filter[, kept_variants], na.rm = T) * 100
  kept_cells <- base::which(mc.c >= mc)
  
  ngt_mutated <- mutated & ngt_filter
  
  ngt_mutated <- ngt_mutated[kept_cells, kept_variants]
  mm.c <- base::colMeans(ngt_mutated, na.rm = T) * 100
  
  mv.f <- (mv.c < mv)
  mm.f <- (mm.c < mm)
  mv.f[!mv.f] <- mm.f
  
  kept_variants <- base::which(!mv.f)
  if (!(length(kept_variants) & length(kept_cells))) {
    stop("All cells/variants are filtered. Try different filtering settings.")
  }
  
  ######## start making changes to tapestri object
  
  filtered_tapestri_object = tapestri_object
  
  if (gt.mask == TRUE) {
    filtered_tapestri_object$variant$features$NGT[!ngt_filter] <- 3
  }
  
  #filter variants from all layers in variants assay
  for(layer in names(filtered_tapestri_object$variant$features)){
    filtered_tapestri_object$variant$features[[layer]] = filtered_tapestri_object$variant$features[[layer]][,kept_variants]
  }
  
  filtered_tapestri_object$variant$annotations = filtered_tapestri_object$variant$annotations[kept_variants,]
  
  #filter cells from all assays
  for(assay in names(filtered_tapestri_object)) {
    for(layer in names(filtered_tapestri_object[[assay]]$features)) {
      
      filtered_tapestri_object[[assay]]$features[[layer]] = filtered_tapestri_object[[assay]]$features[[layer]][kept_cells,]    
    }
  }
  

  str(filtered_tapestri_object, max.level=3, vec.len = 2)
  return(filtered_tapestri_object)
}








#' Filter raw genotypes based on quality
#'
#' @param loom Loom file
#' @param gqc Genotype quality cutoff (default 30)
#' @param dpc Read depth cutoff (default 10)
#' @param afc Allele frequency cutoff (default 20)
#' @param mv Remove variants with < mv of known values (default 50)
#' @param mc Remove variants with < mc of known values (default 50)
#' @param mm Remove variants mutated in < mm of cells (default 1)
#' @param gt.mask mask low quality GT as missing (if GQ/DP/AF lower than cutoff, default FALSE)
#' @return Filtered genotypes
#' @export
filter_variants2 <- function(variant_assay, gqc = 30, dpc = 10, afc = 20, mv = 50, mc = 50, mm = 1, gt.mask = FALSE) {
  
  # variant_assay = variants
  # gqc = 30
  # dpc = 10
  # afc = 20
  # mv = 50
  # mc = 50
  # mm = 1
  # gt.mask = FALSE
  # 
  
  if(variant_assay@assay_name != ASSAY_NAME_VARIANT) {
    stop("Assay is not of name variant.")
  }
  
  data = variant_assay@data_layers
  
  gt <- as.matrix(data$NGT)
  mask <- (gt < 3)
  mutated <- (gt == 1 | gt == 2)
  
  gq <- (data$GQ >= gqc)
  
  
  dp <- data$DP
  ad <- data$AD
  af <- matrix(100, nrow = nrow(dp), ncol = ncol(dp))
  
  af[mutated] <- ad[mutated] * 100 / dp[mutated]
  af[is.na(af)] <- 0
  af <- (af >= afc)
  dp <- (dp >= dpc)
  ngt_filter <- gq & dp & af & mask
  mv.c <- base::colMeans(ngt_filter, na.rm = T) * 100
  kept_variants <- base::which(mv.c >= mv)
  mc.c <- base::rowMeans(ngt_filter[, kept_variants], na.rm = T) * 100
  kept_cells <- base::which(mc.c >= mc)
  
  ngt_mutated <- mutated & ngt_filter
  
  ngt_mutated <- ngt_mutated[kept_cells, kept_variants]
  mm.c <- base::colMeans(ngt_mutated, na.rm = T) * 100
  
  mv.f <- (mv.c < mv)
  mm.f <- (mm.c < mm)
  mv.f[!mv.f] <- mm.f
  
  kept_variants <- base::which(!mv.f)
  if (!(length(kept_variants) & length(kept_cells))) {
    stop("All cells/variants are filtered. Try different filtering settings.")
  }
  
  ######## start making changes to tapestri object
  if (gt.mask == TRUE) {
    variant_assay@data_layers$NGT[!ngt_filter] <- 3
  }
  
  filtered_variant_assay = create_assay(assay_name = ASSAY_NAME_VARIANT, 
                                        cell_annotations = variant_assay@cell_annotations[kept_cells,],
                                        feature_annotations = variant_assay@feature_annotations[kept_variants,]
                                      )
  
  #filter variants from all layers in variants assay
  for(layer in names(variant_assay@data_layers)){
    filtered_variant_assay = add_data_layer(filtered_variant_assay,
                                            layer_name = layer, 
                                            data = variant_assay@data_layers[[layer]][kept_cells,kept_variants])
  }
  
  return(filtered_variant_assay)
}


































