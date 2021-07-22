#' @importFrom BiocFileCache BiocFileCache bfcadd
downloadATAC <- function(df, row, column, bfc){
    if (df[row, column] != ''){
        filename <- bfcadd(bfc, "TestWeb", fpath=df[row,column])
        return(readRDS(filename))
    } else {
        return(NULL)
    }

}

fetchATAC <- function(df, row, sparse){
    #download the data into dataframes
    cache_path <- tempfile()
    bfc <- BiocFileCache(cache_path, ask = FALSE)
    if (sparse == FALSE){
        if (df[row,'dense_matrix_link']==""){
            stop(paste(df[row,'Accession'],
            "has no dense matrix, use sparse=TRUE to download it."))
        }
        expression <- downloadATAC(df, row, 'dense_matrix_link', bfc)
    } else if (sparse == TRUE){
        if (df[row,'sparse_matrix_link']==""){
            stop(paste(df[row,'Accession'],
            "has no sparse matrix, use sparse=FALSE to download it."))
        }
        expression <- downloadATAC(df, row, 'sparse_matrix_link', bfc)
    }
    labels <- downloadATAC(df, row, 'cell_annotation_link', bfc)
    if (!is.null(labels) && length(labels$cell)!=length(colnames(expression))){
        col.num <- which(colnames(expression) %in% labels$cell)
        expression <- expression[,col.num]
    }

    dataset_data_meta <- list(#signatures = sigs,
                        pmid = df[row, 'PMID'],
                        author = df[row, 'Author'],
                        technology = df[row, 'Sequencing_Technology'],
                        score_type = df[row, 'Score_Type'],
                        organism  = df[row, 'Organism'],
                        genome_build = df[row, 'Genome_Build'],
                        cell_categories = df[row,'Broad_Cell_Categories_Present'],
                        tissue_type = df[row,'Tissue_Cell_Type'],
                        disease = df[row,'Disease'],
                        summary = df[row,'Data_Summary'],
                        cells = colnames(expression),
                        matrix_name = df[row, 'Matrix_Names'],
                        #identifiers for peak matrix
                        regions = row.names(expression),
                        accession = df[row, 'Accession'])
    if (is.null(labels)){
        dataset <- SingleCellExperiment(list(counts = expression),
                                            metadata = dataset_data_meta)
    }else{
        dataset <- SingleCellExperiment(list(counts = expression),
                                    colData = data.frame(label=labels[,c("cluster", "cell_label")]),
                                    metadata = dataset_data_meta)
    }


    return(dataset)

}
