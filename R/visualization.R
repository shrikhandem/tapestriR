

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
    clusters = clusters,
    features
  )
  
  ### Make into long format for plotting
  df_long = data_to_plot %>% 
    pivot_longer(-c(clusters), 
                 names_to = 'feature',values_to = 'value') 
  
  
  ### plot
  p = ggplot(data=df_long)
  p = p + geom_violin(aes(x=clusters, y=value, fill=clusters))
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
#' @import tidyverse
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



#' Title
#'
#' @param normalized_reads read counts per amplicon per cell
#' @param clusters cluster labels 
#'
#' @return
#' @export
#' @import tidyverse
#' 
tapestri_ploidy_plot <- function(normalized_reads, clusters) {
  
  data_to_plot = tibble(
    normalized_reads,
    clusters = clusters
  )
  
  df_long = data_to_plot %>% 
    pivot_longer(-c(contains('clusters')), 
                 names_to = 'feature',values_to = 'value') %>% 
    mutate(feature = as_factor(feature)) %>% group_by(feature) %>% mutate(median_count = median(value,na.rm = TRUE)) 
  
  # group by clusters
  p = ggplot(data=df_long, aes(x=feature, y=value))
  
  p  = p + geom_dotplot(aes(fill = clusters),   # Use fill = Species here not in ggplot()
                        binaxis = "y",         # which axis to bin along
                        binwidth = 0.1,        # Minimal difference considered different
                        stackdir = "center",    # Centered
                        alpha=.04
  )
  
  p  = p + stat_summary(fun = median, geom = "crossbar", width = 0.8) + facet_wrap(~clusters, ncol=1) + 
    ylim(c(0,ceiling(max(df_long$median_count))+1))
  p = p + geom_hline(aes(yintercept=2), color="red", size=1, alpha=.5)
  p = p + xlab('') + ylab('')
  p = p + theme_bw()  + theme(legend.position = "none",
                              axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(p)  
  
}



