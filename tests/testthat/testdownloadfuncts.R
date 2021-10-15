
library(scATAC.Explorer)

# should return a list with metadata dataframe inside
test_that("Test metadata table retrieval", {
    expect_type(queryATAC(metadata_only = TRUE), "list")
    expect_s3_class(queryATAC(metadata_only = T)[[1]], "data.frame")
})

# should be able to download small dataset
test_that("Test ability to download dataset", {
    # should retrieve a single dataset
    expect_equal(length(queryATAC(accession = "GSE89362")), 1)
    # downloaded dataset should be a SingleCellExperiment object
    expect_s4_class(queryATAC(accession = "GSE89362")[[1]],
     "SingleCellExperiment")
})