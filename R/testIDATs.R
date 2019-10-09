#' Make sure annotations and so forth are in place before doing a big run...
#'
#' @param samps Two or more rows of the `samps` data.frame with `Basename`
#' @param ... Other arguments to be passed to readMethArrays
#'
#' @return Boolean of if sesamize finished successfully
#'
#' @examples
#'
#' @export
#'
testIDATs <- function(samps,
                      ...) {
    rgSet <- readMethArrays(base = ".", targets = samps[1:2, ], ...)
    colnames(rgSet) <- rgSet$subject
    res <- sesamize(rgSet)
    is(res, "GenomicRatioSet")
}
