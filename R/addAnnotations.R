#' Add annotations to object returned from sesamizeGEO
#'
#' GEOmetadb requires a large memory overhead. This function uses the phenoData
#' data frame to supply the annotation information, rather than using GEOmetadb
#' like addCharacteristics
#'
#' @param x Object returned from sesamizeGEO
#' @param pheno phenoData data frame from GEOquery
#'
#' @return Object returned from sesamizeGEO with added annotations
#'
#' @import minfi
#'
#' @examples
#'
#' @export
#'
addAnnotations <- function(x,
                           pheno) {

    if (is(x, "GenomicRatioSet")) {
        if (!c("subject") %in% names(colData(x))) {
            stop("You need a column named subject in your colData to run this")
        }
        GSMs <- x$subject
    } else {
        GSMs <- x
    }

    if (length(GSMs) != nrow(pheno)) {
        pheno <- subset(pheno, rownames(pheno) %in% GSMs)
    }

    toParse <- pheno[, grep("*:ch1", names(pheno)), drop = FALSE]
    rownames(toParse) <- GSMs
    colnames(toParse) <- elts(colnames(toParse), sep=":", elt=1)
    toParse[] <- lapply(toParse, as.factor)
    toParse$title <- as.character(pheno$title)
    toParse$gsm <- GSMs

    if (is(x, "GenomicRatioSet")) {
        colnames(x) <- toParse$title
        colnames(metadata(x)$SNPs) <- colnames(x)
        for (i in names(toParse)) colData(x)[, i] <- toParse[, i]
        for (i in 1:length(colData(x))) {
            if (class(colData(x)[[i]]) == "factor") {
                names(colData(x)[[i]]) <- colData(x)$gsm
            }
        }
        return(x)
    } else {
        return(toParse)
    }
}
