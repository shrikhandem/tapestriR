

#' Filter raw genotypes based on quality
#' 
#' @param variant_assay Assay object of type dna
#' @param gqc Genotype quality cutoff (default 30)
#' @param dpc Read depth cutoff (default 10)
#' @param afc Allele frequency cutoff (default 20)
#' @param mv Remove variants with < mv of known values (default 50)
#' @param mc Remove variants with < mc of known values (default 50)
#' @param mm Remove variants mutated in < mm of cells (default 1)
#' @param gt.mask mask low quality GT as missing (if GQ/DP/AF lower than cutoff, default FALSE)
#'
#' @return Filtered genotypes
#' @export
filter_variants <- function(variant_assay, gqc = 30, dpc = 10, afc = 20, mv = 50, mc = 50, mm = 1, gt.mask = FALSE) {
  
  # variant_assay = variants
  # gqc = 30
  # dpc = 10
  # afc = 20
  # mv = 50
  # mc = 50
  # mm = 1
  # gt.mask = FALSE

  needed_layers = c("AD","DP","GQ","NGT")
  check_assay = needed_layers %in% names(variant_assay@data_layers)
  if(sum(check_assay)!=4) {
    stop("Assay must containing")
  }
  
  data = variant_assay@data_layers
  
  gt <- as.matrix(data$NGT)
  mask <- (gt < 3)
  mutated <- (gt == 1 | gt == 2)
  
  gq <- (data$GQ >= gqc)
  
  
  dp <- as.matrix(data$DP)
  ad <- as.matrix(data$AD)
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


































