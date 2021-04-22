# scATAC-Explorer

## SKELETON README, STILL NEED TO EXPLAIN MULTIPLE MATRICES, AND ADD IMAGES  

## Introduction

scATAC-Explorer is a curated collection of publicly available scATAC-seq (Single Cell Assay for Transposase-Accessible Chromatin using sequencing) datasets. It aims to provide a single point of entry for users looking to study chromatin accessibility at the single-cell level. 

Users can quickly search available datasets using the metadata table, and then download the datasets they are interested in for analysis. Optionally, users can save the datasets for use in applications other than R. 

This package will improve the ease of studying and integrating scATAC-seq datasets. Developers may use this package to obtain data for analysis of multiple tissues, diseases, cell types, or developmental stages. It can also be used to obtain data for validation of new algorithms. 


## Installation
``` 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
library(devtools)
install_github("shooshtarilab/scATACdb")
```

# Tutorial

## Exploring available datasets

Start by exploring the available datasets through metadata.

```
> res = queryATAC(metadata_only = TRUE)
```

This will return a list containing a single dataframe of metadata for all available datasets. View the metadata with `View(res[[1]])` and then check `?queryATAC` for a description of searchable fields.

Note: in order to keep the function's interface consistent, `queryATAC` always returns a list of objects, even if there is only one object. You may prefer running `res = queryATAC(metadata_only = TRUE)[[1]]` in order to save the dataframe directly.

![Screenshot of the metadata table](docs/metadata.png)

The `metatadata_only` argument can be applied alongside any other argument in order to examine only datasets that have certain qualities. You can, for instance, view only carcinoma datasets by using 

```
> res = queryATAC(disease = 'Carcinoma', metadata_only = TRUE)[[1]]
```

![Screenshot of the metadata table](docs/bc_metadata.png)

| Search Parameter | Description                                     | Examples                |
| ---------------- | ----------------------------------------------- | ----------------------- |
| geo_accession    | Search by GEO accession number                  | GSE129785, GSE89362     |
| has_cell_types   | Filter by presence of cell-type annotations     | TRUE, FALSE             |
| has_clusters     | Filter by presence of cluster results           | TRUE, FALSE             |
| disease          | Search by disease                               | Carcinoma, Leukemia     |
| author           | Search by first author                          | Satpathy, Cusanovich    |
| journal          | Search by publication journal                   | Science, Nature, Cell   |
| year             | Search by year of publication                   | <2015, >2015, 2013-2015 |
| pmid             | Search by PubMed ID                             | 27526324, 32494068      |
| sequence_tech    | Search by sequencing technology                 | 10x Genomics Chromium   |
| organism         | Search by source organism                       | Mus musculus            |
| sparse           | Return expression in sparse matrices            | TRUE, FALSE             |

#### Searching by year

In order to search by single years and a range of years, the package looks for specific patterns. '2013-2015' will search for datasets published between 2013 and 2015, inclusive. '<2015' will search for datasets published before or in 2015. '>2015' will search for datasets published in or after 2015.


### Getting your first dataset

Once you've found a field to search on, you can get your data. 

```
> res = queryATAC(geo_accession = "GSE131688")
```

This will return a list containing dataset GSE131688. The dataset is stored as a `SingleCellExperiment` object, with the following metadata list:

#### Metadata
| Attribute     | Description |
| ------------- | --------------------------------------------------------------- |
| cells         | A list of cells included in the study |
| regions       | A list of genomic regions (peaks) included in the study |
| pmid          | The PubMed ID of the study |
| technology    | The sequencing technology used |
| genome_build  | The genome build used for data generation |
| score_type    | The type of scoring or normalization used on the counts data |
| organism      | The type of organism from which cells were sequenced |
| author        | The first author of the paper presenting the data |
| disease       | The diseases sampled cells were sampled from |
| summary       | A broad summary of the study conditions the sample was assayed from |
| geo_accession | The GEO accession ID for the dataset |

#### Accessing data

To access the expression data for a result, use
```
> View(counts(res[[1]]))
```
![Screenshot of the metadata table](docs/GSE72056_expression.png)

Cell type labels are stored under `colData(res[[1]])` for datasets for which cell type labels are available.

To access metadata for a dataset, use
```
> metadata(res[[1]])
```
Specific metadata entries can be accessed by specifying the attribute name, for instance
```
> metadata(res[[1]])$pmid
```


### Example: Returning all datasets with cell-type labels

Say you want to measure the performance of cell-type classification methods. To do this, you need datasets that have the true cell-types available. 
```
> res = queryATAC(has_cell_types = TRUE)
```
This will return a list of all datasets that have true cell-types available. You can see the cell types for the first dataset using the following command:
```
> View(colData(res[[1]]))
```
![Screenshot of the cell type labels](docs/GSE72056_labels.png)

The first column of this dataframe contains the cell barcode or cell ID, the second contains the cell type, and the third contains the cluster assignment if available. 

## Saving Data

To facilitate the use of any or all datasets outside of R, you can use `saveATAC()`. `saveATAC` takes two parameters, one a `ATAC_data` object to be saved, and the other the directory you would like data to be saved in. Note that the output directory should not already exist.

To save the data from the earlier example to disk, use the following commands.

```
> res = queryATAC(geo_accession = "GSE131688")[[1]]
> saveATAC(res, '~/Downloads/GSE131688')
[1] "Done! Check ~/Downloads/GSE131688 for files"
```
The result is three files (a counts .mtx file, a peak region .tsv file, and a cell ID/Barcodes .tsv file) that can be used in other programs. In the future we will support saving in other formats.

#### NOTE: need to confirm this, might be fine for scATAC-seq data as .mtx file is already a sparse format, so no conversion to dense neccessary.
NOTE: `saveATAC` is currently not compatible with sparse datasets. This is due to the size of some datasets and the memory required to convert them to a dense matrix that can be written to a csv file. To save the elements of a sparse object, use `write.table()` and `as.matrix(counts(res))`, keeping in mind that doing this with some of the larger datasets may cause R to crash.

![Screenshot of the saveTME files](docs/saveTME_files.png)


## System Requirements

While many of the datasets included in this package are small enough to be loaded and stored, even as dense matrices, on machines with an 'average' amount of memory (8-32gb), there are a few larger datasets that cannot be fully manipulated without a significant amount of memory. With this in mind, we recommend using `sparse = TRUE` when possible and using a system with at least 64gb of RAM for full functionality.

If you are experience crashes due to memory limitations, try using `sparse = TRUE` or grabbing datasets individually using the `geo_accession` parameter.
