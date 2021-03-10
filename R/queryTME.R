#' A function to query TME datasets available in this package
#'
#' This function allows you to search and subset included TME datasets. 
#' A list of tme_data objects matching the provided options will be returned, 
#' if queryTME is called without any options it will retrieve all available datasets. 
#' This should only be done on machines with a large amount of ram (>64gb) because some datasets are quite large.
#' In most cases it is recommended to instead filter databases with some criteria.
#' @param geo_accession Search by geo accession number. Good for returning individual datasets
#' @param author Search by the author who published the dataset
#' @param journal Search by the journal the dataset was published in.
#' @param year Search by exact year or year ranges with '<', '>', or '-'. For example, you can return datasets newer than 2013 with '>2013'
#' @param pmid Search by Pubmed ID associated with the study. Good for returning individual datasets
#' @param sequence_tech Search by sequencing technology used to sample the cells.
#' @param score_type Search by type of score (TPM, FPKM, raw count)
#' @param has_clusters Return only those datasets that have clustering results available, or only those without (TRUE/FALSE)
#' @param has_truth Return only those datasets that have cell-type annotations available, or only those without annotations
#' @param organism Search by source organism used in the study, for example human or mouse.
#' #TODO update docs
#' @param genome_build
#' @param category
#' @param tissue
#' @param disease
#' @param metadata_only Return rows of metadata instead of actual datasets. Useful for exploring what data is available without actually downloading data. Defaults to FALSE
#' @param sparse Return expression as a sparse matrix. 
#'                  Uses less memory but is less convenient to view, recommended only if encounter memory issues with dense data. Defaults to FALSE.
#' @keywords tumour
#' @importFrom methods new
#' @importFrom Matrix Matrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @return A list containing a table of metadata or 
#' one or more SingleCellExperiment objects
#'
#' @examples
#' 
#' ## Retrieve the metadata table to see what data is available
#' res <- queryTME(metadata_only = TRUE)
#' 
#' ## Retrieve a filtered metadata table that only shows datasets with 
#' ## cell type annotations and cell type gene signatures
#' res <- queryTME(has_truth = TRUE, has_signatures = TRUE, metadata_only = TRUE)
#' 
#' ## Retrieve a single dataset identified from the table
#' res <- queryTME(geo_accession = "GSE72056")

queryTME <- function(geo_accession=NULL,
                    author=NULL, 
                    journal=NULL, 
                    year=NULL, 
                    pmid=NULL, 
                    sequence_tech=NULL, 
                    score_type=NULL,
                    has_clusters=NULL,
                    has_truth=NULL,
                    organism=NULL,
                    genome_build=NULL, #TODO add this code
                    category=NULL, #TODO add code
                    tissue=NULL, #TODO add code
                    disease=NULL, #TODO add code

                    metadata_only=FALSE,
                    sparse = FALSE){
    df <- tme_meta
    if (!is.null(geo_accession)) {
        df <- df[df$accession == geo_accession,]
    }
    if (!is.null(author)) {
        df <- df[toupper(df$author) == toupper(author),]
    }
    if (!is.null(journal)) {
        df <- df[toupper(df$journal) == toupper(journal),]
    }
    if (!is.null(year)) {
        year <- gsub(' ', '', year)
        #check greater than
        if (gregexpr('<', year)[[1]][[1]] == 5 || gregexpr('>',year)[[1]][[1]]==1){
            year <- sub('>','',year)
            year <- sub('<','',year)
            df <- df[df$year>=year,]
        
        #check between
        }else if (grepl('-',year,fixed=TRUE)){
            year <- strsplit(year,'-')[[1]]
            df <- df[df$year>=year[[1]]&df$year<=year[[2]],]
        
        #check less than
        }else if (gregexpr('>', year)[[1]][[1]] == 5 || gregexpr('<',year)[[1]][[1]]==1){
            year <- sub('>','',year)
            year <- sub('<','',year)
            df <- df[df$year<=year,]
        
        #check equals
        }else{
            df <- df[df$year==year,]
        }
    }
    if (!is.null(pmid)) {
        df <- df[df$PMID == pmid,]
    }
    if (!is.null(sequence_tech)) {
        #TODO
        df <- df[toupper(df$sequencing_tech) == toupper(sequence_tech),]
    }
    if (!is.null(score_type)) {
        #TODO eventually this will become a way to select which type of score you want to 
        # download since we will store multiple types
        df <- df[toupper(df$score_type) == toupper(score_type) ,]
    }
    if (!is.null(has_clusters)) {
        if (has_clusters) {
            df <- df[df$clusters == 'Y' ,]
        }else if (!has_clusters) {
            df <- df[df$clusters == 'N' ,]
        }
    }
    if (!is.null(has_truth)) {
        if (has_truth) {
            df <- df[df$cell_labels == 'Y', ]
        }else if (!has_truth) {
            df <- df[df$cell_labels == 'N', ]
        }
    }
    if (!is.null(organism)) {
        df <- df[toupper(df$Organism) == toupper(organism),]
    }
    if (!is.null(genome_build)) {
        df <- df[toupper(df$Genome_Build) == toupper(genome_build),]
    }
    if (!is.null(category)) {
        #TODO this column is a list of strings separated by "+", 
        # need to actually read it as a list
        df <- df[toupper(df$broad_cell_categories) == toupper(category),]
    }
    if (!is.null(tissue)) {
        #TODO also a list
        df <- df[toupper(df$tissue_type) == toupper(tissue),]
    }
    if (!is.null(disease)) {
        #TODO also a list
        df <- df[toupper(df$Disease) == toupper(disease),]
    }

    if (metadata_only) {
        df[,c('signature_link', 'expression_link', 'truth_label_link','sparse_expression_link')] <- list(NULL)
        return(list(df))
    } else {
        df_list <- list()
        for (row in seq_len(nrow(df))){
            df_list[[row]] <- fetchTME(df, row, sparse)
        }
        return(df_list)
    }

    return(list(df))

}