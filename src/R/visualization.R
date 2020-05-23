
######## 2d projection Theme
projection_theme <- theme_bw() + theme(axis.line=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.text.y=element_blank(),
                                       axis.ticks=element_blank(),
                                       axis.title.x=element_blank(),
                                       axis.title.y=element_blank(),
                                       panel.background=element_blank(),
                                       panel.grid.major=element_blank(),
                                       panel.grid.minor=element_blank())

x_y_theme <- theme_bw() 



#' Scatter plot wrapper for ggplots that makes plotting Tapestri object easier
#'
#' @param data Tapestri object
#' @param x plot of x
#' @param y plot of y
#' @param color_by what features of color by 
#' @param colorby expression or genotypes values for a feature
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' p <- multiassay_scatterplot(model, "umap", cluster, 'Test umap plot')
#' }
scatterplot <- function(data, x, y, color_by) {
  
  
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


violinplot <- function(data, clusters , features) {
  
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
  p = p + x_y_theme + theme(legend.position = "none",
                            axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(p)
  
}


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
