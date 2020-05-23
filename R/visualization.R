

#' Simple theme to show axis-less projection of umaps
#'
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
#' @param collapse_zygosity 
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
