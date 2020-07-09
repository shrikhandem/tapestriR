# Overview

Using the Mission Bio Tapestri Platform, reagents, and the `tapestriR` package, we explore single nucleotide variants (SNVs), protein expression, and copy number variations (CNVs), including loss of heterozygosity (LOH), all of which play a large role in cancer evolution and contribute to cancer heterogeneity. More details can be found in the [Application Note] (https://missionbio.com/protein-reva-application-note-form/). 

The `tapestriR` package is an exploratory data analysis tool. Refer to the [README.md] (https://github.com/MissionBio/tapestriR/blob/master/README.md) file for instructions on how to install the `tapestriR` package. The latest version can be found at https://github.com/MissionBio/tapestriR.


In the vignettes, we show a few example methods to normalize, cluster, and visualize multi-omics data. Users can also explore different methods available in the R community to further explore and interpret their data. For an easy-to-use data analysis and visualization solution, check out [Tapestri Insights] (https://portal.missionbio.com/).

# Loading Data

The main value of this package is to load data from these sources and create an R object for use in additional analysis and visualization:
* Tapestri H5 file - enables analysis of DNA variants, CNVs, and protein expression.
* Loom file - enables analysis of DNA variants only.
* Exported data from Tapestri Insights - enables analysis of DNA variants only.

# Objects
The Multiomics Object is a structured way to store all the raw and normalized data as well as the analysis performed on the data. The object ensures that all the multi-omics data is properly aligned for cells and features (variants, protein counts, etc.) of the cells.

* row = cells
* columns = features

## Multiomics Object

A Multiomics Object stores a group of assays performed on 1 or more samples. 

```
<root>
+-- metadata                          Experiment metadata. Required to have experiment_name.
+-- assays
          +-- dna_variants            DNA Assay Object
          +-- dna_read_counts         DNA Read Counts Assay Object
          +-- protein_read_counts     Protein Read Counts Assay Object
          +--
+-- cell_annotations                  Annotation of each cell. Required to have barcode and sample name columns.

```


## Assay Object

An Assay Object stores raw and processed data for 1 or more samples. 

* An assay can contain multiple data layers and analysis layers.
* Data layers must have the same cell barcodes and features.
* Analysis layers must have the same cell barcodes but can have different features.
** Typically, store clustering and dimension reduction values in the analysis layers. 

### `dna_variants` Assay Object

```
<root>
+-- assay_name                        Name of the assay. Same as one of the layers in the assays slot of the Multiomics Object.
+-- data_layers                       1 or more layers of data. Each layer must have the same dimensions and matching cell barcodes (rows) and features (columns names). 
                +--NGT                Genotype layer, 0=WT, 1=HET, 2=HOM, 3=unknown
                +--AD
                +--GQ
                ...                   Additional layers can be added. Must have the same cell barcodes and features. 
+-- analysis_layers                   Additional layers after the data is processed. Must have the same number of cells (rows). 
                    +--               Examples: umap projection of each cell or kmeans cluster labels.
+-- feature_annotations               Annotations of features. Must have an id column.
+-- cell_annotations                  Annotations of cells. Required to have barcode and sample name columns. 
```

### `dna_read_counts` Assay Object

```
<root>
+-- assay_name                        Name of the assay. Same as one of the layers in the assays slot of the Multiomics Object.
+-- metadata                          Assay metadata
+-- data_layers                       1 or more layers of data. Each layer must have the same dimensions and matching cell barcodes (rows) and features (columns names). 
                +--read_counts        Raw read counts per DNA
                +--normalized         Normalized read counts
+-- analysis_layers                   Additional layers after the data is processed. Must have the same number of cells (rows). 
                    +--
+-- feature_annotations               Annotations of features. Must have an id column.
+-- cell_annotations                  Annotations of cells. Required to have barcode and sample name columns. 
```

### `protein_read_counts` Assay Object

```
<root>
+-- assay_name                        Name of the assay. Same as one of the layers in the assays slot of the Multiomics Object.
+-- metadata                          Assay metadata
+-- data_layers                       1 or more layers of data. Each layer must have the same dimensions and matching cell barcodes (rows) and features (columns names). 
                +--read_counts        Raw read counts per protein
                +--normalized         Normalized read counts
+-- analysis_layers                   Additional layers after the data is processed. Must have the same number of cells (rows). 
                    +--
+-- feature_annotations               Annotations of features. Must have an id column.
+-- cell_annotations                  Annotations of cells. Required to have barcode and sample name columns. 
```


