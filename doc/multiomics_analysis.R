## ----include=FALSE------------------------------------------------------------

# title: "Multi-omics Analysis"
# output: rmarkdown::html_vignette
# vignette: >
#   %\VignetteIndexEntry{Multi-omics Analysis}
#   %\VignetteEngine{knitr::rmarkdown}
#   \usepackage[utf8]{inputenc}
# editor_options: 
#   chunk_output_type: console


library(TapestriR)
require(rhdf5)
library(tidyverse)
library(ggplot2)
library(plotly)
require(patchwork)

# umap library
library(uwot)

#normalization 
library(SciViews)

# clustering
require(scran)


ASSAY_NAME_VARIANT = 'dna'
ASSAY_NAME_READ_COUNT = 'cnv'
ASSAY_NAME_PROTEIN = 'protein'



## ----load_data----------------------------------------------------------------

filename <- system.file("extdata", "ABseq021.h5", package = "TapestriR")

# ideally would just start by loading filtered H5, but for now will send filtering list to function to mimic
variants = read_h5(filename = filename, assay_name = ASSAY_NAME_VARIANT, min_mutation_rate = 0.01)

filtered_variants = filter_variants(variants)

vaf=round(filtered_variants@data_layers$AD/filtered_variants@data_layers$DP, 3)
vaf[is.na(vaf)] <- 0
filtered_variants = add_data_layer(filtered_variants,'VAF',vaf)

protein = read_h5(filename,assay_name = ASSAY_NAME_PROTEIN)
protein_counts_norm = protein@data_layers$read_counts %>% clr_by_feature() %>% as_tibble(rownames = NA)
protein = add_data_layer(protein,'normalized',protein_counts_norm)

cnv = read_h5(filename,assay_name = ASSAY_NAME_READ_COUNT)

ab21 = create_moo(experiment_name = basename(filename), cell_annotations = filtered_variants@cell_annotations)
ab21 = add_assay(ab21,filtered_variants)
ab21 = add_assay(ab21,cnv,keep_common_cells = TRUE)
ab21 = add_assay(ab21,protein, keep_common_cells = TRUE)

experiment = ab21

## ----fig.width = 12, fig.height = 8-------------------------------------------


##################
# select the Proteins to plot on X and Y
##################

#protein_x = 'CD34'
#protein_y = 'CD38'

protein_x = 1
protein_y = 2

##################
# select 1 or more features to color by
# color_by should be a vector of column header you want to color by
##################

#all proteins
color_by = experiment@protein@data_layers$normalized

#select a few proteins
#color_by =  experiment@protein@data_layers$normalized %>% select('CD110','CD117')

#select a few variant
#color_by =  experiment@dna@data_layers$NGT %>% select(1:10) %>% mutate_all(as_factor) %>% mutate_all(recode_genotypes)


p  = tapestri_scatterplot(x = experiment@protein@data_layers$normalized[[protein_x]], 
                 y= experiment@protein@data_layers$normalized[[protein_y]], 
                 color_by = color_by)
p = p + xlab(protein_x) + ylab(protein_y)
p = p + theme_bw()  + scale_colour_gradient2(low="yellow", mid='grey', high="darkred") 
p

#select a few variant
color_by =  experiment@dna@data_layers$NGT %>% select(1:10) %>% mutate_all(as_factor) %>% mutate_all(recode_genotypes)


p  = tapestri_scatterplot(x = experiment@protein@data_layers$normalized[[protein_x]], 
                 y= experiment@protein@data_layers$normalized[[protein_y]], 
                 color_by = color_by)
p = p + xlab(protein_x) + ylab(protein_y)
p = p + theme_bw()  
p



## -----------------------------------------------------------------------------


## ----cluster snv, cache=TRUE--------------------------------------------------

##################
# do demension reduction and clustering a few different ways
##################

#dimensional reduction using umap
set.seed(111)
umap_values <- umap(experiment@dna@data_layers$VAF, scale=TRUE, metric="manhattan", init="laplacian", pca=20) 

umap_layer = tibble(    
      x = umap_values[,1],
      y = umap_values[,2]
)

experiment@dna = add_analysis_layer(assay = experiment@dna, layer_name = 'umap_vaf', umap_layer)

######
# Hold all the different customer labels in a single structure
#####
cluster_by = experiment@dna@analysis_layers$umap_vaf
clusters = list()

#### do the clustering

for(i in 2:5) {
  kmean_values <- kmeans(cluster_by, i ,iter.max=500)
  clusters[[paste0('umap.kmean.cluster.',i)]] = as_factor(kmean_values$cluster)

  kmean_values <- kmeans(cluster_by, i ,iter.max=500)
  clusters[[paste0('feature.kmean.cluster.',i)]] = as_factor(kmean_values$cluster)

}

graph_values <- buildSNNGraph(t(cluster_by), k=150)
louvain_clust <- igraph::cluster_louvain(graph_values)$membership

clusters[['umap.louvain.cluster']] = as_factor(louvain_clust)

graph_values <- buildSNNGraph(t(cluster_by), k=150)
louvain_clust <- igraph::cluster_louvain(graph_values)$membership

clusters[['features.louvain.cluster']] = as_factor(louvain_clust)


#############
## Add cluster labels to analysis data structure
#############
experiment@dna = add_analysis_layer(assay = experiment@dna, layer_name = 'umap_vaf_clusters', as_tibble(clusters))



## ----fig.width = 12, fig.height = 8-------------------------------------------

p  = tapestri_scatterplot(
                 x = experiment@dna@analysis_layers$umap_vaf$x, 
                 y= experiment@dna@analysis_layers$umap_vaf$y, 
                 color_by = experiment@dna@analysis_layers$umap_vaf_clusters)
p = p + umap_theme() + ggtitle('umap_vaf_clusters')
p


## ----fig.width = 12, fig.height = 12------------------------------------------
 
 #%>% select(!contains('chr2:198267'))

color_by = experiment@dna@data_layers$NGT %>% select(1:2) %>% mutate_all(as_factor) %>% mutate_all(recode_genotypes)


p  = tapestri_scatterplot(
                 x = experiment@dna@analysis_layers$umap_vaf$x, 
                 y= experiment@dna@analysis_layers$umap_vaf$y, 
                 color_by = color_by)
p = p + umap_theme() 
p = p + ggtitle('Color by Genotypes')
ggplotly(p)



## ----fig.width = 12, fig.height = 8-------------------------------------------


color_by = experiment@dna@analysis_layers$umap_vaf_clusters$umap.kmean.cluster.2

p  = tapestri_scatterplot(
                 x = experiment@dna@analysis_layers$umap_vaf$x, 
                 y= experiment@dna@analysis_layers$umap_vaf$y, 
                 color_by = color_by)
p = p + umap_theme() + ggtitle('umap_vaf_clusters') + theme(legend.position = 'none')


v = tapestri_violinplot(clusters = color_by,
               features = experiment@protein@data_layers$normalized)
v = v + theme_bw() + theme(legend.position = "none",
                            axis.text.x = element_text(angle = 90, hjust = 1))
  
## pathwork magic
p / v



## ----fig.width = 12, fig.height = 8-------------------------------------------

p  = tapestri_scatterplot(
                 x = experiment@dna@analysis_layers$umap_vaf$x, 
                 y= experiment@dna@analysis_layers$umap_vaf$y, 
                 color_by = experiment@protein@data_layers$normalized)
p = p + umap_theme() + scale_colour_gradient2(low="yellow", mid='grey', high="darkred") 
p = p + ggtitle('Color by Proteins')
p


## ----cache=TRUE---------------------------------------------------------------

##################
# do demension reduction and clustering a few different ways
##################

#dimensional reduction using umap
set.seed(111)
umap_values <- umap(experiment@protein@data_layers$normalized, scale=TRUE, metric="manhattan", init="laplacian", pca=20) 

umap_layer = tibble(    
      x = umap_values[,1],
      y = umap_values[,2]
)

experiment@protein = add_analysis_layer(assay = experiment@protein, layer_name = 'umap', umap_layer)
  
######
# Hold all the different customer labels in a single structure
#####
cluster_by = experiment@protein@analysis_layers$umap
clusters = list()

#### do the clustering

for(i in 2:5) {
  kmean_values <- kmeans(cluster_by, i ,iter.max=500)
  clusters[[paste0('umap.kmean.cluster.',i)]] = as_factor(kmean_values$cluster)

  kmean_values <- kmeans(cluster_by, i ,iter.max=500)
  clusters[[paste0('feature.kmean.cluster.',i)]] = as_factor(kmean_values$cluster)

}

graph_values <- buildSNNGraph(t(cluster_by), k=150)
louvain_clust <- igraph::cluster_louvain(graph_values)$membership

clusters[['umap.louvain.cluster']] = as_factor(louvain_clust)

graph_values <- buildSNNGraph(t(cluster_by), k=150)
louvain_clust <- igraph::cluster_louvain(graph_values)$membership

clusters[['features.louvain.cluster']] = as_factor(louvain_clust)


#############
## Add cluster labels to analysis data structure
#############
experiment@protein = add_analysis_layer(assay = experiment@protein, layer_name = 'clusters', as_tibble(clusters))



## ----fig.width = 12, fig.height = 8-------------------------------------------


p  = tapestri_scatterplot(
                 x = experiment@protein@analysis_layers$umap$x, 
                 y= experiment@protein@analysis_layers$umap$y, 
                 color_by = experiment@protein@analysis_layers$clusters)

p = p + xlab('') + ylab('')
p = p + umap_theme() 
p


## ----fig.width = 12, fig.height = 8-------------------------------------------

p_umap  = tapestri_scatterplot(
                 x = experiment@protein@analysis_layers$umap$x, 
                 y= experiment@protein@analysis_layers$umap$y, 
                 color_by = experiment@protein@analysis_layers$clusters$features.louvain.cluster) + 
  xlab('') + ylab('') + umap_theme() 

p_violin = tapestri_violinplot(
           clusters = experiment@protein@analysis_layers$clusters$features.louvain.cluster,
           features = experiment@protein@data_layers$normalized)

p_umap / p_violin


## ----fig.width = 12, fig.height = 8-------------------------------------------

p  = tapestri_scatterplot(
                 x = experiment@protein@analysis_layers$umap$x, 
                 y= experiment@protein@analysis_layers$umap$y, 
                 color_by = experiment@protein@data_layers$normalized)
p = p + umap_theme() + scale_colour_gradient2(low="yellow", mid='grey', high="darkred") 
p


## ----fig.width = 12, fig.height = 8, eval=FALSE-------------------------------
#  
#  data_to_plot = tibble(
#    experiment@dna@analysis_layers$clusters,
#    analysis_df$cnv$features$read_counts
#  )
#  
#  df_long = data_to_plot %>%
#    pivot_longer(-c(contains('cluster')),
#                 names_to = 'feature',values_to = 'value') %>%
#    mutate(feature = as_factor(feature))
#  
#  # group by clusters
#  p = ggplot(data=df_long
#             #%>% filter(feature %in% colnames(analysis_df$amplicon$features)[30:60])
#       )
#  p = p + geom_boxplot(aes(x=feature, y=value, fill=umap.kmean.cluster.2), outlier.shape = NA, notch = TRUE )
#  #p = p + facet_wrap(~umap.kmean.cluster, ncol = 1)
#  #p = p + geom_smooth(aes(x=as.numeric(feature), y=value, color=umap.kmean.cluster), method = 'loess',se=FALSE)
#  p = p + xlab('') + ylab('')
#  p = p + theme_bw()  + theme(legend.position = "none",
#                axis.text.x = element_text(angle = 90, hjust = 1))
#  p
#  

## -----------------------------------------------------------------------------


