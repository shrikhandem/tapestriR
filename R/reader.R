# global variables
if (getRversion() >= "2.15.1") utils::globalVariables(c("variant", "alt", "start", "end"))

# class union
setClassUnion(name="MultiAssayOrNull", c("MultiAssay", "NULL"))

#' Read NGT (genotypes) data
#'
#' @param ngt_file NGT file exported from Insights
#' @return MultiSample genotypes data
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' ngt <- read_genotypes("NGT.csv")
#' }
read_genotypes <- function(ngt_file) {
  type = "snv"
  ngt_data <- read.table(ngt_file, sep=",", header=T, check.names=F)
  filtered <- TRUE
  normalized <- TRUE
  samples <- split(ngt_data, ngt_data$Sample)
  dna_analytes <- list()
  for (i in c(1:length(samples))) {
    dna <- samples[[i]]
    name <- as.character(dna$Sample[1])
    data <- dplyr::select(dna, -dplyr::one_of("Sample", "Cell", "Subclone"))
    col_names <- gsub(":", "\\.", colnames(data))
    col_names <- gsub("/", "\\.", col_names)
    colnames(data) <- col_names
    rownames(data) <- dna$Cell
    barcodes <- dna$Cell
    clones <- NULL
    if ('Subclone' %in% colnames(dna)) {
      clones <- dplyr::select(dna, c("Sample", "Cell", "Subclone"))
    } else {
      clones <- dplyr::select(dna, c("Sample", "Cell"))
      clones$Subclone <- clones$Sample
    }
    analyte <- methods::new(Class = 'Analyte',
      type = type,
      data = data,
      name = name,
      barcodes = barcodes,
      features = colnames(data),
      filtered = filtered,
      normalized = normalized,
      clones = clones)
    dna_analytes[[name]] <- analyte
  }
  #ngt <- ngt_data[, c(4:ncol(ngt_data))]
  #clones <- ngt_data[,c(1,3,2)]
  #clones[,3] <- as.character(clones[,3])
  #clones[stringr::str_detect(clones[,3], "clones"),3] <- "Unknown"
  #identity <- ngt_data[,1]
  #multisample <- methods::new(Class = "MultiSample", samples=unique(identity), 
  #  type="snv", data = ngt, clones=clones, features = colnames(ngt),
  #  barcodes = ngt_data[,3], identity = identity)
  return(dna_analytes)
}


#' Filter protein analytes data to remove invalid barcodes
#'
#' @param analyte Protein analyte
#' @param barcodes valid barcodes
#' @return Analyte Filtered protein analyte
#' @export
#' @examples
#' \dontrun{
#' analyte <- filter_protein_barcodes("NGT.csv")
#' }
filter_protein_barcodes <- function(analyte, barcodes) {
  if (analyte$type != "prot") {
    stop("Please specify protein analyte")
  }
  protein <- analyte$counts
  if (!any(stringr::str_detect(barcodes, "-\\d$"))) {
    rownames(protein) <- gsub("-.*", "", rownames(protein))
  }

  prot_barcodes <- rownames(protein)
  comm_barcodes <- barcodes[!is.na(match(barcodes, prot_barcodes))]
  if(length(comm_barcodes) == 0) {
    stop("No matching barcodes found")
  }
  protein <- protein[comm_barcodes, ]
  analyte$barcodes <- comm_barcodes
  analyte$counts <- protein
  analyte$filtered <- TRUE

  return(analyte)
}


#' Create protein analyte object
#'
#' @param name Sample name
#' @param filename input file
#' @return Analyte Analyte object
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' analyte <- read_protein_counts("input.tsv", "sample1")
#' }
read_protein_counts <- function(filename, name=NULL) {
  type = "prot"
  if (is.null(name)) {
    name <- unlist(strsplit(basename(filename), "\\."))[1]
  }
  data <- read.table(filename, header=T, sep="\t")
  data <- as.data.frame(reshape2::acast(data, cell_barcode ~ ab_description, value.var="raw", fill=0))
  col_names <- sort(colnames(data))
  data <- data[, col_names]

  filtered <- NULL
  normalized <- FALSE
  if (type == "prot") {
    filtered <- FALSE
  }
  barcodes = rownames(data)
  clones <- c()
  clones$Cell <- barcodes
  clones$Sample <- name
  clones$Subclone <- name
  clones <- as.data.frame(do.call(cbind, clones))
  analyte <- methods::new(Class = 'Analyte', type = type, counts = data, name = name,
    barcodes = rownames(data),features = colnames(data), 
    filtered = filtered, normalized = normalized, clones = clones)
  return(analyte)
}

# Create MultiSample object
#
# @param analytes Analytes of same type from different samples
# @param type Type of analyte
# @return MultiSample object
# @importFrom methods new
# @export
# @examples
# \dontrun{
# multi_prot <- create_multi_samples(analytes, "prot")
# }
create_multi_samples <- function(analytes, type){

  raw_data <- data.frame()
  norm_data <- data.frame()
  clones <- data.frame()
  features <- c()
  
  for (i in c(1:length(analytes))) {
    analyte <- suppressWarnings(filter_slots_by_type(analytes[[i]], type))

    # skip if analyte not found
    if (is.null(analyte)) { next }

    if (analyte$type != type) {
      stop("Please specify analytes of same type")
    }

    if (is.null(analyte$data) && type != "snv") {
      stop("Please run normalization on analytes")
    }

    if (!is.null(analyte$filtered) && analyte$filtered == FALSE && type == "prot") {
      stop("Please filter protein data")
    }
    if (!is.null(analyte$normalized) && analyte$filtered == FALSE && type %in% c("dna", "prot")) {
      stop("Please normalize dna/protein data")
    }

    data <- analyte$counts
    if (type == "snv") {
      data <- analyte$data
    }

    if (is.null(features)) {
      features <- colnames(data)
    }

    features <- intersect(features, colnames(data))
    data <- data[, features]
    data <- tibble::rownames_to_column(data, "Cell")
    data$Samples <- analyte$name
    if (nrow(raw_data) != 0) {
      raw_data <- raw_data[, c("Cell", features, "Samples")]
      norm_data <- norm_data[, features]
    }
    raw_data <- rbind(raw_data, data)
    clones <- rbind(clones, analyte$clones)
    if (type == "snv") {
      norm_data <- rbind(norm_data, analyte$data[, features])
    } else {
      if (nrow(analyte$data) == 0) {
        warning("No normalized data found, using raw values")
        norm_data <- rbind(norm_data, analyte$counts[, features])
      } else {
        norm_data <- rbind(norm_data, analyte$data[, features])
      }
    }
  }

  if (nrow(raw_data) == 0) {
    #warning(sprintf("No analytes for %s assay are found.", type))
    return(NULL)
  }

  multisample <- methods::new(Class = "MultiSample", samples=unique(raw_data$Samples), 
    type=type, data = norm_data, barcodes = raw_data$Cell, clones=clones,
    features = colnames(norm_data), identity = raw_data$Samples)

  return(multisample)
}


#' Create MultiSample object
#'
#' @param name name of the analysis
#' @param analytes Analytes of same type from different samples
#' @return TapestriR object
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' multi_prot <- merge_samples(analytes)
#' }
merge_samples <- function(name, analytes){
  multisample_objects <- list()
  for (assay in c('dna', 'prot', 'snv')) {
    object <- create_multi_samples(analytes, assay)
    if (!is.null(object)) {
      multisample_objects[[assay]] <- object
    }
  }

  dna <- multisample_objects[['dna']]
  snv <- multisample_objects[['snv']]
  prot <- multisample_objects[['prot']]
  
  version <- packageVersion(pkg = 'tapestri')

  multi_assay <- methods::new(
    Class = "TapestriR",
    name = name,
    dna=dna,
    snv=snv,
    prot=prot,
    version=version)

  return(multi_assay)
}

#' Filter object by type
#'
#' @param analytes List of MultiSample objects
#' @param type Type of object to be filtered
#' @return MultiSample object
#' @export
#' @examples
#' \dontrun{
#' multi_prot <- filter_objects_by_type(analytes, "prot")
#' }
filter_objects_by_type <- function(analytes, type) {

  for (analyte in analytes) {
    if (analyte$type == type) {
      return(analyte)
    }
  }

}


#' Create multiassay-multisample object
#'
#' @param name Name of the analysis
#' @param analytes MultiSample objects from various assays
#' @return TapestriR object
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' tapestrir <- merge_assay("project1", analytes)
#' }
merge_assays <- function(name, analytes) {

  dna = filter_objects_by_type(analytes, 'dna')
  snv = filter_objects_by_type(analytes, 'snv')
  prot = filter_objects_by_type(analytes, 'prot')
  if(is.null(snv) && is.null(prot) && is.null(dna)) {
    stop("Please provide atleast one valid analyte")
  }

  barcodes <- NULL
  dna_samples <- NULL
  snv_samples <- NULL
  prot_samples <- NULL
  dna_barcodes <- NULL
  snv_barcodes <- NULL
  prot_barcodes <- NULL

  sname <- NULL

  if (!is.null(dna)) {
    dna_barcodes <- dna$barcodes
    sname <- dna$name
    dna_samples <- dna$name
    dna_barcodes <- gsub("-.*", "", dna_barcodes)
  } 
  if (!is.null(snv)) {
    sname <- snv$name
    snv_barcodes <- snv$barcodes
    snv_samples <- snv$name
    snv_barcodes <- gsub("-.*", "", snv_barcodes)
  } 
  if (!is.null(prot)) {
    sname <- prot$name
    prot_barcodes <- prot$barcodes
    prot_samples <- prot$name
    prot_barcodes <- gsub("-.*", "", prot_barcodes)
  } 
  if (is.null(dna_barcodes) && is.null(snv_barcodes) && is.null(prot_barcodes)) {
    stop("No barcodes found. Exiting...")
  }

  common <- Reduce(intersect, Filter(Negate(is.null),
    list(dna_samples, snv_samples, prot_samples)))

  barcodes <- Reduce(intersect, Filter(Negate(is.null),
    list(dna_barcodes, snv_barcodes, prot_barcodes)))

  if (length(barcodes) == 0) {
    stop("No matching barcodes found. Exiting...")
  }

  filtered_samples <- list()
  for (analyte in c(dna, snv, prot)) {
    if (analyte$type != "snv" && analyte$normalized == FALSE) {
      stop("Please run normalization on counts data(dna/protein)")
    }
    data_barcodes <- rownames(analyte$data)
    data_barcodes <- gsub("-.*", "", data_barcodes)
    rownames(analyte$data) <- data_barcodes

    if (!is.null(analyte$data)) {
      analyte$data <- analyte$data[barcodes,]
    }
    if (!is.null(analyte$counts) && nrow(analyte$counts) != 0) {
      rownames(analyte$counts) <- data_barcodes
      analyte$counts <- analyte$counts[barcodes,]
    }
    if (!is.null(analyte$clones) && nrow(analyte$clones) != 0) {
      clones <- analyte$clones
      clones$Cell <- gsub("-.*", "", clones$Cell)
      clones <- clones[clones$Cell %in% barcodes, ]
      analyte$clones <- clones
    }
    filtered_samples[[analyte$type]] <- analyte
  }

  if (length(common) == 0) {
    stop("Sample ID's do not match. Please use correct Sample ID's")
  }

  multi_assay <- methods::new(
    Class = "MultiAssay",
    name=sname,
    dna=filtered_samples[['dna']],
    snv=filtered_samples[['snv']],
    prot=filtered_samples[['prot']],
    barcodes=barcodes)

  return(multi_assay)
}

#' Create dna(read counts) analyte object
#'
#' @param name Sample name
#' @param files input files
#' @return Analyte Analyte object
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' analyte <- read_dna_counts("sample", c("tube0.tsv", "tube1.tsv"))
#' }
read_dna_counts <- function(name, files) {
  if(is.null(name)) {
    stop("Please specify sample name")
  }
  data <- read_barcodes(files)
  barcodes <- rownames(data)
  clones <- c()
  clones$Cell <- barcodes
  clones$Sample <- name
  clones$Subclone <- name
  clones <- as.data.frame(do.call(cbind, clones))
  analyte <- methods::new(Class = 'Analyte', type="dna", counts=data, name=name,
   features=colnames(data), barcodes = barcodes, filtered = TRUE,
   normalized = FALSE, clones = clones)
  return(analyte)
}


#' Convert dataframe to Analyte object
#'
#' @param data data
#' @param type assay type
#' @param name sample name
#' @return Analyte Analyte object
#' @importFrom methods new
#' @export
#' @examples
#' \dontrun{
#' genotypes <- extract_genotypes(loom, barcodes, gt.mask=TRUE)
#' analyte <- convert_to_analyte(genotypes, "snv", "sample")
#' }
convert_to_analyte <- function(data, type, name) {
  # add gene names as suffix to variant ids
  if (type == "snv") {
    col_names <- colnames(data)
    col_names <- annotate_variants_with_genes(col_names)
    colnames(data) <- col_names
  }
  data <- as.data.frame(data)

  analyte <- methods::new(Class = 'Analyte',
    type = type,
    data = data,
    name = name,
    barcodes = rownames(data),
    features = colnames(data),
    filtered = TRUE,
    normalized = TRUE,
    clones = data.frame())
  return(analyte)
}


#
# Annotate variant ids with genes and add gene name as prefix to variant id
# @param col_names variant ids
# @return variants annotated col names
#
annotate_variants_with_genes <- function(col_names) {
  out <- as.data.frame(col_names)
  colnames(out) <- "variant"
  out$name <- out$variant
  out$variant <- gsub('\\.\\.\\.', '\\.N\\.', out$variant)
  out <- out %>% tidyr::separate(variant, into=c("chr", "start", "ref", "alt")) %>% 
    dplyr::mutate(alt = stringr::str_length(alt)) %>% 
    dplyr::mutate(end=as.numeric(start) + as.numeric(alt) - 1) %>% 
    dplyr::select(chr, start, end, name)
  out <- annotate_with_genes(out, col_names)
  out$id[!is.na(out$annot.symbol)] <- 
    paste(out$annot.symbol[!is.na(out$annot.symbol)], 
      out$id[!is.na(out$annot.symbol)], sep=".")

  return(out$id)
}




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
read_loom <- function(filename, min_mutation_rate=0.05) {
  
  # filename <- "data/PE11.cells.loom"
  # min_mutation_rate = 0.1
  # experiment_name = basename(filename)
  
  h5f = rhdf5::H5Fopen(filename)
  
  features=as_tibble(h5f$row_attrs)
  cells = h5f$col_attrs
  
  
  layers = h5f$layers
  ngt = h5f$matrix
  layers$NGT = ngt
  cell_annotations = tibble(cell_id=1:nrow(ngt))
  
  
  mutation_rate = apply(ngt,2, FUN = function(x) {
    sum(x %in% c(1,2)) / ncol(ngt)
  })
  mutated_variants = mutation_rate > min_mutation_rate
  filtered_features = features[mutated_variants,]
  
  assay_name = ASSAY_NAME_VARIANT
  
  assay = create_assay(assay_name = assay_name,
                       cell_annotations = cell_annotations,
                       feature_annotations = filtered_features)
  
  for(layer in names(layers)) {
    data = layers[[layer]][,mutated_variants]
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay,layer_name = layer, data = data)
  }
  
  
  rhdf5::H5Fclose(h5f)
  return(assay)
}


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
h5_reader <- function(filename, assay_name, min_mutation_rate = 0.01) {
  
  # filename <- "~/Google Drive/launches/r_package/insights_v3/data/ABseq021.h5"
  # assay_name=ASSAY_NAME_VARIANT
  # layer = 'NGT'
  # min_mutation_rate = 0.1

  h5f = rhdf5::H5Fopen(filename)
  
  if(assay_name == ASSAY_NAME_VARIANT) {
    
    # filter data to make it managable in R
    ngt = h5f$assays$dna$layers$NGT
    mutation_rate = apply(ngt,1, FUN = function(x) {
      sum(x %in% c(1,2)) / ncol(ngt)
    })
    filters = mutation_rate > min_mutation_rate
  } else {
    filters = TRUE
  }
  
  h5_assay = h5f$assays[[assay_name]]
  
  filtered_features <- as_tibble(h5_assay$ca)[filters,]
  cell_annotations <- as_tibble(h5_assay$ra) %>% mutate_all(as.character)
  
  assay = create_assay(assay_name = assay_name,
                     cell_annotations = cell_annotations,
                     feature_annotations = filtered_features)

  
  for(layer in names(h5_assay$layers)) {
    data = t(h5_assay$layers[[layer]][filters, ])
    colnames(data) <- filtered_features$id
    assay = add_data_layer(assay = assay, layer_name = layer, data = data)
  }

  rhdf5::H5Fclose(h5f)
  devnull <- base::gc()
  return(assay)     
}
