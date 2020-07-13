#' Define subclones based on NGT values
#' 
#' Similar to how Tapestri Insights v2.2 calculates subclones 
#'
#' @param assay Assay Object for DNA variants
#' @param small_subclone threshold (%) for considering subclones too small
#' @param ignore_zygosity collapse HET/HOM to just mutated
#'
#' @return
#' @export
#'
define_subclones <- function(assay, small_subclone = 0.01, ignore_zygosity=FALSE){
  
  ngt = experiment$assays$dna$data_layers$NGT

    if(ignore_zygosity) {
    ngt[ngt ==2] = 1
  }
  
  df = ngt %>% group_by_all() %>% mutate(
    cells = n(),
    subclone = paste('subclone',group_indices())
  ) %>% ungroup()
  
  #df %>% arrange(-cells) %>% select(contains('subclone'),cells) 
  
  df <- df %>% mutate(
    subclone_label = case_when (
      apply(df == 3, 1, any) ~ 'missing NGT subclone',
      (cells/nrow(experiment$assays$dna)) < small_subclone ~ 'small subclone',
      TRUE ~ subclone
    ))
  
  # df %>% group_by(subclone_label) %>% summarise(cells=n()) %>% arrange(-cells)
  return(df)
}
