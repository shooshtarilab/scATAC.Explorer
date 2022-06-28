rm(list=ls())
library(Matrix)
library(tibble)
library(tidyr)

# defining paths
matrix_dir = "../Downloaded/"
barcode.path <- "../Processed/GSE96769_scATACseq_barcodes.txt"
features.path <- paste0(matrix_dir, "Peaks/GSE96769_PeakFile_20160207.bed")
matrix.path <- paste0(matrix_dir, "Cell_Matrix/GSE96769_scATACseq_counts.txt")

# reading in matrix, peak, barcode files
mat = read.table(matrix.path, 
                 header = FALSE,
                 sep="\t",
                 stringsAsFactors = FALSE)

feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)[,1:3]

barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

# convert into our sparse matrix format from csr matrix
mat = sparseMatrix(
  i = mat$V1,
  j = mat$V2,
  x = mat$V3,
  dims = c(max(mat$V1),max(mat$V2))
)
# might not need dim declaration, ncol and nrow line up with peaks and barcodes

# check dimensions
if(ncol(mat) != nrow(barcode.names)){print("ERROR: Number of matrix columns is not equal to the number of barcode names")}
if(nrow(mat) != nrow(feature.names)){print("ERROR: Number of matrix rows is not equal to the number of feature names")}

# column names of matrix are cell barcodes
colnames(mat) = barcode.names[,1]

# merge peak info columns together, using "-" as a seperator
feature.names <- data.frame(feature.names) %>% unite("regions", 1:3 , sep = '-')

# add in regions as rowname for matrix
feature.names <- as.vector(as.matrix(feature.names))
rownames(mat) <- feature.names

mat <- Matrix(mat, sparse = TRUE)

# save matrix
output.path = "../Processed/GSE96769_scATACseq_Matrix_Sparse.rds"

saveRDS(mat, file = output.path)



