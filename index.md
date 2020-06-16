# Overview

Using the Mission Bio Tapestri Platform, reagents and the *tapestriR* package, we explore Single Nucleotide Variants (SNVs), protein expression and Copy Number Variations (CNVs), including loss of heterozygosity (LOH), which plays a large role in cancer evolution and contributes to cancer heterogeneity. More details can be found in the [Application Note] (https://missionbio.com/protein-reva-application-note-form/). 

The *tapestri* R package is an exploratory data analysis tool. Refer to the README.md file for instructions on how to install the *tapestri* R package. 

Here we show a few example methods to normalize, cluster and visualize multi-omics data. It is left to the user to explore different methods avaliable in the R community to better explore and interpret their data. 

# Objects

The main value of this package is to load data from these sources and create an R object for use in additional analysis and visualiation.
* Tapestri H5 file
* loom file
* Exported data from Tapestri Insights


The multiomics object is a structured way to store all of the raw and normalized data as well as the analysis performed on these data. The object ensures that all of the multiomics data is properly aligned for cells and features of the cells.

* row = cells
* columns = features


### Multiomics Object
```
<root>
+-- metadata                          Experiment meta data. Required to have experiment_name
+-- assays
          +-- dna_variants            DNA Assay object
          +-- dna_read_counts         DNA read counts Assay object
          +-- protein_read_counts     Protein Read counts Assay object
          +--
+-- cell_annotations                  Annotation of each cell, Required to have barcode and sample name columns

```

Assays must have identical cell barcodes. Each assay will have different featuers. 

#### `dna_variants` Assay Object
```
<root>
+-- assay_name                        Name of Assay. Same as one of the layers in the assays slot of Multiomics object
+-- data_layers                       1 or more layers of data. each layer must have same dimensions and matching cell barcodes (rows) and features (columns names) 
                +--NGT                Genotype layer. 0=WT, 1=HET, 2=HOM, 3=unknown
                +--AD
                +--GQ
+-- analysis_layers                   Additional layers after data is processed. Must have same number of cells (rows) 
                    +--               Examples: umap projection of each cell, or kmeans cluster labels
+-- feature_annotations               Annotations of features. Must have id column that is same as column names in data layers.
+-- cell_annotations                  Annotations of cells. Required to have barcode and sample name columns 
```

* An assay can multiple data layers and analysis layers
* data layers must have same cell barcodes and features
* analysis layers must have same cell barcodes, but can have different features
** Typically store clustering and dimension reduction values in the analysis layers 

#### `protein_read_counts` Assay Object
```
<root>
+-- assay_name                        Name of Assay. Same as one of the layers in the assays slot of Multiomics object
+-- metadata                          Name of Assay. Same as one of the layers in the assays slot of Multiomics object
+-- data_layers                       1 or more layers of data. each layer must have same dimensions and matching cell barcodes (rows) and features (columns names) 
                +--read_counts        raw read counts per protein
                +--normalized         normalized read counts
                +--
+-- analysis_layers                   Additional layers after data is processed. Must have same number of cells (rows) 
                    +--
+-- feature_annotations               Annotations of features. Must have id column that is same as column names in data layers.
+-- cell_annotations                  Annotations of cells. Required to have barcode and sample name columns 
```
