#' Retrieve phenoData data frame from GEO series
#'
#' Given a GEO series accessed by GEOquery::getGEO(), return the phenoData data
#' frame
#'
#' @param series GEO series accessed by GEOquery::getGEO()
#' @param gse String with the GEO series name
#' @param gpl String with the GEO platform name
#'
#' @return phenoData data frame
#'
#' @examples
#'
#' @export
#'
getPhenoData <- function(series,
                         gse,
                         gpl) {
    
    if (substr(gse, start = 1, stop = 3) != "GSE") {
        gse <- paste0("GSE", gse)
    }
    if (substr(gpl, start = 1, stop = 3) != "GPL") {
        gpl <- paste0("GPL", gpl)
    }

    test_name <- paste0(gse,"-",gpl,"_series_matrix.txt.gz", sep = "")

    if (test_name %in% names(series)) {
        pheno <- pData(series[[test_name]])
    } else {
        test_name <- paste0(gse,"_series_matrix.txt.gz", sep = "")
        pheno <- pData(series[[test_name]])
    }

    message("phenoData taken from: ", test_name)
    return(pheno)
}
