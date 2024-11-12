
#' A function to save a scATAC-seq dataset stored in a SingleCellExperiment
#'
#' This function allows you to save the counts,
#' peaks, cell ID's/barcodes, and any cell clustering data to disk in csv format. It takes two options:
#' an object to save and a directory to save in. Multiple files will be created in
#' the provided output directory, one for each type of data available in the scATAC_data object
#' (counts, cell ID/Barcode, peak regions, cell type/cluster annotations).
#' @param object The SingleCellExperiment object to be written to disk, this should be an individual dataset returned by queryATAC.
#' @param outdir The directory to save the data in, the directory should not exist yet.
#' @keywords scATAC-seq
#' @importFrom methods is
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment colData
#' @importFrom Matrix Matrix
#' @importFrom utils write.csv
#' @importFrom zellkonverter writeH5AD
#' @export
#' @return Nothing
#'
#' @examples
#'
#' # Retrieve a previously identified dataset (see queryATAC) and save it to disk
#' res <- queryATAC(accession = 'GSE89362')[[1]]
#' \dontshow{
#'          #res <- SingleCellExperiment(list(counts=matrix()))
#'          tdir = tempdir()
#'          output_directory_name = file.path(tdir, 'save_tme_data')}
#' saveATAC(res, output_directory_name)
#'
saveATAC <- function(object, outdir, format = "mtx") { 
    if (!is(object, "SingleCellExperiment")) {
        stop('object parameter must be of type SingleCellExperiment')
    }
    if (file.exists(outdir)) {
        stop('outdir must not be an existing directory')
    } else {
        dir.create(outdir)
    }
    if (format == "mtx") {
        expr_name <- file.path(outdir,
                                paste(object@metadata$accession,
                                "matrix.mtx",
                                sep = '_'))
        cellID_name <- file.path(outdir,
                                paste(object@metadata$accession,
                                "barcodes.tsv",
                                sep = '_'))
        peaks_name <- file.path(outdir,
                                paste(object@metadata$accession,
                                "peaks.tsv",
                                sep = '_'))
        label_name <- file.path(outdir,
                                paste(object@metadata$accession,
                                "cell_types_and_clusters.csv",
                                sep = '_'))

        # will always have matrix, cellID, and peaks
        Matrix::writeMM(SingleCellExperiment::counts(object), file = expr_name)
        write(colnames(object), file = cellID_name)
        write(row.names(object), file = peaks_name)
        # only write cluster/celltype data if we have it
        if (length(colnames(colData(object)) > 0)) {
            write.csv(colData(object), file = label_name)
        }
    }
    if (format == "h5ad") {
        writeH5AD(object, "h5ad_file.h5ad")

    }
    print(paste('Done! Check', outdir, 'for files', sep = ' '))
}
