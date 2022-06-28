rm(list=ls())
library(Matrix)
library(tibble)
library(tidyr)

setwd("C:/Users/AG/OneDrive/Year 4/Thesis/Processing Scripts/")

# defining paths to downloaded files in market matrix format
# downloaded files are in format of .mtx file containing read counts (0 or 1) and two text files 
# containing cell ID's and peak region names
matrix_dir <- "Example Data/"
barcode.path <- paste0(matrix_dir, "GSE126074_P0_BrainCortex_SNAREseq_chromatin.barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "GSE126074_P0_BrainCortex_SNAREseq_chromatin.peaks.tsv.gz")
matrix.path <- paste0(matrix_dir, "GSE126074_P0_BrainCortex_SNAREseq_chromatin.counts.mtx.gz")

# reading in matrix, peak, barcode files
mat <- readMM(file = matrix.path)
feature.names <- read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names <- read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

# checking to make sure the length of the counts matrix is the same length as the number of cells in the cell ID
# and that the number of rows in the counts matrix matches the number of peak regions 
if(ncol(mat) != nrow(barcode.names)){print("ERROR: Number of matrix columns is not equal to the number of barcode names")}
if(nrow(mat) != nrow(feature.names)){print("ERROR: Number of matrix rows is not equal to the number of feature names")}

# column names of matrix are cell barcodes 
colnames(mat) <- barcode.names$V1

# want to format peak regions as "chromosme-start-end"
# replace underscores in original peak regions with -
feature.names[,1] <- gsub(":", "-", feature.names[,1])

# assign rownames of counts matrix to be the peak regions
feature.names <- as.vector(as.matrix(feature.names))
rownames(mat) <- feature.names

# make sure matrix is in sparse format at end
mat <- Matrix(mat, sparse = TRUE)

# make sure sparse matrix is in dgCMatrix format, somtimes previous code will convert to dgTMatrix, which we don't want
mat <- as(mat, "dgCMatrix")

# save matrix
output.path = "/Example Data/GSE126074_P0_BrainCortex_SNAREseq_chromatin_Matrix.rds"
saveRDS(mat, file = output.path)

