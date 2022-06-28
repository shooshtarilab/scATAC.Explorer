rm(list=ls())
library(Matrix)
library(tidyr)

# defining paths to already processed sparse matrix
input.path = "../Processed/GSE129785_scATAC-TME-All-Matrix-Sparse.rds"
output_folder = "../Processed/"

# reading in processed sparse matrix
mat = readRDS(input.path)

# separate the peak region info into seperate chr, start, end columns
feature.names = as.data.frame(rownames(mat))
feature.names = separate(feature.names,1, c("Chr","Start","End"), sep = "-")[,1:3]
feature.names = cbind(feature.names,Region = rownames(mat))

# remove rownames from sparse matrix and convert to dense
rownames(mat) = NULL
mat = as.matrix(mat)
# add the feature.names column to mat
mat = cbind(feature.names,mat)

# checking to make sure the number of peak regions is the same as the number of rows in the dense matrix
if(ncol(mat) != ncol(mat)){print("ERROR: Number of sparse matrix columns is not equal to the number of cells in original matrix")}
if(nrow(mat) != nrow(feature.names)){print("ERROR: Number of matrix rows is not equal to the number of feature names")}

# defining paths and saving
dense_matrix.path <- paste0(output_folder, "GSE129785_scATAC-TME-All-Matrix-Dense.rds")
saveRDS(mat, file = dense_matrix.path)
