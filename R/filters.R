
#' Read Tapestri multi-omics h5 file
#'
#' @param filename h5 file
#' @param min_mutation_rate only read variants with mutations rate higher then treshold. Primarly done to reduce size of data in memory
#' @return Tapestri multi-omics object
#' @export
#' @examples
#' \dontrun{
#' tapestri_raw = h5_reader(filename,min_mutation_rate = 0.1)
#' }
h5_reader <- function(filename, min_mutation_rate = 0.01) {
  
  # filename <- "~/Google Drive/launches/r_package/insights_v3/data/ABseq021.h5"
  # assay_name='dna'
  # feature_filter = colnames(genotypes.zyg)
  # layer = 'NGT'
  # min_mutation_rate = 0.1
  
  h5f = rhdf5::H5Fopen(filename)
  
  # filter data to make it managable in R
  ngt = h5f$assays$dna$layers$NGT
  mutation_rate = apply(ngt,1, FUN = function(x) {
    sum(x %in% c(1,2)) / ncol(ngt)
  })
  mutated_variants = mutation_rate > min_mutation_rate
  
  barcodes <- unlist(h5f$assays$dna$ra)
  
  tapestri_object = tibble(cell_id = barcodes)
  
  ###### special load of DNA to reduce size
  assay = h5f$assays[['dna']]
  
  col_names <- unlist(assay$ca$id)[mutated_variants]
  
  
  layers = tibble(cell_id = barcodes)
  for(layer in names(assay$layers)) {
    
    data = t(assay$layers[[layer]][mutated_variants, ])
    colnames(data) <- col_names
    layers[[layer]] = as_tibble(data)
  }
  
  tapestri_object[['variant']] = tibble(features = layers[,-1])
  
  
  ###### load rest of assays
  assay_names <- c("cnv", "protein")
  
  for(assay_name in assay_names) {
    
    assay = h5f$assays[[assay_name]]
    
    col_names <- unlist(assay$ca$id)
    barcodes <- unlist(assay$ra)
    
    layers = tibble(cell_id = barcodes)
    for(layer in names(assay$layers)) {
      
      data = t(assay$layers[[layer]])
      colnames(data) <- col_names
      layers[[layer]] = as_tibble(data)
    }
    tapestri_object[[assay_name]] = tibble(features = layers[,-1])
  }
  
  # map_dfc(tapestri_object, unlist)
  str(tapestri_object, max.level=6, vec.len = 2)
  # str(tapestri_object, max.level=3, vec.len = 2)
  
  # tapestri_object = as_tibble(tapestri_object)
  rhdf5::H5Fclose(h5f)
  devnull <- base::gc()
  return(tapestri_object)     
}


#' Add features to tapestri object
#'
#' @param tapestri_object tapestri object to update
#' @param assay_name feature belongs to assay variants, cnv, proteins 
#' @param feature_name name of feature
#' @param new_data new date to add to taesptri object
#'
#' @return updated tapestri object
#' @export
#'
add_feature <- function(tapestri_object, assay_name, feature_name, new_data) {
  
  # tapestri_object = analysis_df
  # new_data = protein_counts_norm
  # assay_name = 'protein'
  # feature_name='normalized'
  # 
  
  if (!assay_name %in% names(tapestri_object)) {
    stop(paste0("Assay ", assay_name, ' does not exist.'))
  }
  
  current_feature_names = colnames(tapestri_object[[assay_name]]$features[[1]])
  new_feature_names = colnames(new_data)
  
  if (!all.equal(current_feature_names, new_feature_names)) {
    stop('New feature does not have same column names.')
  }
  
  
  tapestri_object[[assay_name]]$features[[feature_name]] = tibble(new_data)
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
  
  tapestri_object = tapestri_raw
  gqc = 30
  dpc = 10
  afc = 20
  mv = 50
  mc = 50
  mm = 1
  gt.mask = FALSE
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
  
  #filter cells from all assays
  filtered_tapestri_object = filtered_tapestri_object[kept_cells,]
  str(filtered_tapestri_object, max.level=6, vec.len = 2)
  return(filtered_tapestri_object)
}































