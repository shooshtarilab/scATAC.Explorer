rm(list=ls())
library(Matrix)
library(tibble)
library(tidyr)
library(rhdf5)
library(Seurat)

# defining paths
matrix_dir = "../Downloaded/"
matrix.path <- paste0(matrix_dir, "Cell_Matrix/GSM4226310_Gfi1-R412X-R412X-filtered_peak_bc_matrix.h5")

# reading in already made sparse matrix from h5 file using Seurat
mat <- Read10X_h5(matrix.path)

# changing from seurat chr:start-end format to chr-start-end format
rownames(mat) <- gsub(":", "-", rownames(mat))

# save matrix
output.path = "../Processed/GSM4226310_Gfi1-R412X-R412X-filtered_peak_bc_matrix_Sparse.rds"

saveRDS(mat, file = output.path)

