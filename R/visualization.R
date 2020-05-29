

#' Simple theme to show axis-less projection of umaps
#' @import ggplot2
#' @return ggplot theme
#' @export
umap_theme <- function() {
  theme_bw() + 
  theme(axis.line=element_blank(),
       axis.text.x=element_blank(),
       axis.text.y=element_blank(),
       axis.ticks=element_blank(),
       axis.title.x=element_blank(),
       axis.title.y=element_blank(),
       panel.background=element_blank(),
       panel.grid.major=element_blank(),
       panel.grid.minor=element_blank()
   )

}


#' Scatter plot wrapper for ggplots that makes plotting Tapestri object easier
#'
#' @param x plot of x
#' @param y plot of y
#' @param color_by what features of color by 
#'
#' @import tidyverse
#' @import ggplot2
#' @return ggplot object
#' @export
tapestri_scatterplot <- function(x, y, color_by) {
  
  
  # data = analysis_df
  # x = analysis_df$protein$features$CD110
  # y = analysis_df$protein$features$CD33
  # 
  # color_by = analysis_df$protein$features
  
  ##################
  # create a new tibble with columns of data you want plot, 
  # make sure the x, y and color_by are in this new tibble
  ##################
  data_to_plot = tibble(
    x = x,
    y = y,
    color_by
  )
  
  ### Make into long format for plotting
  df_long = data_to_plot %>% 
    pivot_longer(-c(x,y), 
                 names_to = 'feature',values_to = 'value') 
  
  
  ### plot
  p = ggplot(data=df_long)
  p = p + geom_point(aes(x = x, y = y, color=value), alpha = 0.5, size=0.8)
  p = p + facet_wrap(~feature)
  return(p)
  
}


#' violin plot wrapper for ggplots that makes plotting Tapestri object easier
#'
#' @param clusters clusters to split data by
#' @param features features to plot as violin graph
#'
#' @import tidyverse
#' @import ggplot2

#' @return
#' @export
tapestri_violinplot <- function(clusters , features) {
  
  # data = analysis_df
  # clusters = analysis_df$protein$clusters$umap.kmean.cluster.2
  # features = analysis_df$protein$features
  
  
  ##################
  # create a new tibble with columns of data you want plot, 
  # make sure the x, y and color_by are in this new tibble
  ##################
  data_to_plot = tibble(
    cluster = clusters,
    features
  )
  
  ### Make into long format for plotting
  df_long = data_to_plot %>% 
    pivot_longer(-c(cluster), 
                 names_to = 'feature',values_to = 'value') 
  
  
  ### plot
  p = ggplot(data=df_long)
  p = p + geom_violin(aes(x=cluster, y=value, fill=cluster))
  p = p + facet_wrap(~feature,nrow =1) + coord_flip()
  p = p + xlab('') + ylab('')
  
  return(p)
  
}


#' Recode NGT values HET, HOM, WT or MUT/WT if zygosity is ignored
#'
#' @param x NGT vector
#' @param collapse_zygosity default TRUE (WT,MUT). if FALSE (WT, HET, HOM)
#'
#' @return
#' @export
#'
recode_genotypes <- function(x,collapse_zygosity=TRUE) {
  if (collapse_zygosity) {
    fct_recode(x,
               'WT' = '0',
               'MUT' = '1',
               'MUT' = '2',
               'unknown' = '3'
               
    )
  } else {
    fct_recode(x,
               'WT' = '0',
               'HET' = '1',
               'HOM' = '2',
               'unknown' = '3'
               
    )
    
  }
}


#' Plot heatmap
#'
#' @param analytes multi assay object
#' @param cluster cluster labels
#' @param snv_features snv features
#' @param prot_features prot features
#' @param groupby group by gene/chr and select variant with highest mutation rate
#' @return heatmap object
#' @export
#' @examples
#' \dontrun{
#' p <- multiassay_heatmap(analytes, cluster, snv_features=c("DNMT3A", "JAK2"),
#'    prot_features=c("CD3", "CD34"), groupby="gene")
#' }
multiassay_heatmap <- function(moo, cluster) {
  clazz <- class(analytes)[1]
  if (!(clazz %in% c("MultiAssay", "TapestriR"))) {
    stop("Please specify multi assay object")
  }
  
  if (!is.null(groupby) && groupby != "gene") {
    stop("Please provide valid groupby option.")
  }
  
  snv <- filter_slots_by_type(analytes, 'snv')
  prot <- filter_slots_by_type(analytes, 'prot')
  
  snv.h <- NULL
  protein.h <- NULL
  
  if (!is.null(snv)) {
    genotypes <- snv$data
    
    #genotype <- genotypes[, snv_features]
    
    genotypes.mat <- genotype_matrix(genotypes, groupby)
    
    #genotypes.mat <- genotypes.mat[order(match(rownames(genotypes.mat), snv@clones$Cell)), ]
    
    legend_params <- list(title = "Genotype", at = c(0, 1, 2, 3),
                          border = "black", labels = c("WT", "HET", "HOM", "Missing"), 
                          title_gp = grid::gpar(fontsize = 10, fontface = "bold"), 
                          labels_gp = grid::gpar(fontsize = 7), legend_height = 5, legend_width = 5, 
                          grid_height = grid::unit(5, "mm"), grid_width = grid::unit(5, "mm"))
    col_fun <- c("grey", "red", "darkred", "white")
    
    if (!(is.null(snv_features))) {
      if (is.null(groupby)) {
        snv_features <- tapestri.cnv.loh.sort_by_loc(snv_features)
      }
      genotypes.mat <- genotypes.mat[, snv_features]
    }
    
    snv.h <- ComplexHeatmap::Heatmap(as.matrix(genotypes.mat), name = "GT",
                                     column_order = colnames(genotypes.mat), split=factor(cluster),
                                     show_row_names=FALSE, row_title_gp = grid::gpar(fontsize = 6),
                                     col=col_fun, heatmap_legend_param=legend_params,
                                     column_names_gp = grid::gpar(fontsize=8),
                                     show_column_dend=FALSE)
  }
  
  if (!is.null(prot)) {
    protein <- prot$data
    if (!is.null(prot_features)) {
      protein <- protein[, prot_features]
    }
    
    protein.h <- ComplexHeatmap::Heatmap(as.matrix(protein), name = "Normalized Protein Counts", 
                                         split=factor(cluster), show_row_names=FALSE,  row_title_gp = grid::gpar(fontsize = 5), 
                                         column_names_gp = grid::gpar(fontsize=8))
  }
  
  return(protein.h + snv.h)
}


