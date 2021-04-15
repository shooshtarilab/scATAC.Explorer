#' @importFrom BiocFileCache BiocFileCache bfcadd
downloadTME <- function(df, row, column, bfc){
    if (df[row, column] != ''){
        filename <- bfcadd(bfc, "TestWeb", fpath=df[row,column])
        return(readRDS(filename))
    } else {
        return(NULL)
    }

}

fetchTME <- function(df, row, sparse){
    #download the data into dataframes
    cache_path <- tempfile()
    bfc <- BiocFileCache(cache_path, ask = FALSE)
    if (sparse == FALSE){
        expression <- downloadTME(df, row, 'dense_matrix_link', bfc)
    } else if (sparse == TRUE){
        expression <- downloadTME(df, row, 'sparse_matrix_link', bfc)
    }
    labels <- downloadTME(df, row, 'cell_annotation_link', bfc)
    if (!is.null(labels) && length(labels$cell)!=length(colnames(expression))){
        col.num <- which(colnames(expression) %in% labels$cell)
        expression <- expression[,col.num]
    }
    #TODO for some reason this is failing, what is different about the metadata
    #sigs <- downloadTME(df, row, 'signature_link', bfc)

    tme_data_meta <- list(#signatures = sigs,
                        pmid = df[row, 'PMID'],
                        author = df[row, 'author'],
                        technology = df[row, 'sequencing_tech'],
                        score_type = df[row, 'score_type'],
                        organism  = df[row, 'Organism'],
                        genome_build = df[row, 'genome_build'],
                        cell_categories = df[row,'broad_cell_categories'],
                        tissue_type = df[row,'tissue_type'],
                        disease = df[row,'disease'],
                        summary = df[row,'Data.Summary']
                        cells = colnames(expression),

                        #TODO maybe figure out how to make this a dataframe with
                        #the first few columns if a dataset has multiple
                        #identifiers for each gene
                        regions = row.names(expression),
                        geo_accession = df[row, 'accession'])
    if (is.null(labels)){
        tme_dataset <- SingleCellExperiment(list(counts = expression),
                                            metadata = tme_data_meta)
    }else{
        tme_dataset <- SingleCellExperiment(list(counts = expression),
                                    colData = data.frame(label=labels$truth),
                                    metadata = tme_data_meta)
    }


    return(tme_dataset)

}
