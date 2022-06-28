rm(list=ls())
library(Matrix)
library(tibble)
library(tidyr)

# defining paths
cluster_dir = "../Downloaded/ClusterData/"

cluster.files.list <- list.files(path = cluster_dir)

cluster.path <- paste0(cluster_dir, file.name)

# reading in cluster info
cluster.info <- read.table(file = cluster.path, sep = "\t", header = FALSE)

#rename columns
colnames(cluster.info) <- c("cell","cluster","cell_label")

# split file name from path and name file [filename]_Cluster.rds
cluster.filename <- strsplit(file.name,"\\.")[[1]][2]
cluster.filename <- paste0(cluster.filename,"_Cluster.rds")

# save dataframe
output.path = paste0("../Processed/",cluster.filename)

saveRDS(cluster.info, file = output.path)




