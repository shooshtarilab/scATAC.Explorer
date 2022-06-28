---
name: New Dataset Submission
about: Suggest a new dataset for inclusion to package
title: "[Dataset] <name>"
labels: addition
assignees: agibsonk

---

Note: R scripts are available under the .github of the package repository for R scripts to process some of the common scATAC-seq data formats for submission to the package. 

**Dataset description**

Please provide a description of the dataset, such as what type of cells were sequenced, the number of cells, etc.

**Link to the reference for the dataset**

Please provide a DOI, article reference, pubmed ID, or other method of linking to the paper.


**Link to the hosted dataset**

This should be a link to the element(s) of the dataset. It can be a link to any file-hosting platform as long as the files are in the format outlined below. 

- cell-by-peak matrix in a R dgCMatrix object saved as a .rds file 
    - The rownames of the matrix should be the genomic regions in "chromosome-start-end" format
    - The column names of the matrix should be the names of the individual cells (this doesn't need to have a particular format, often they are barcodes)
    - The contents of the matrix should be the chromatin accessibility data
- Cell-type labels should be R dataframes saved in .rds files
    - Rownames of the dataframe should be the cell name
    - The first column should be named "cluster" and contain cluster assignments for the specified cell in the same row
    - The second column should be named "cell_label" and contain the cell type label for the cell in the same row

**Metadata Table**

Please fill all columns of the new_scATAC_metadata.csv and attached the completed document here.