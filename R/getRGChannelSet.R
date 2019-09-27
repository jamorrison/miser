#' Assemble a suitable RGChannelSet for sesamizing
#' 
#' Note: this should probably go away in favor of signalSet parallel processing
#' 
#' @param subjects The names of each subject, or NULL to autodetect
#'                  (DEFAULT: NULL)
#' @param frags Which elements to extract for the array Basename (DEFAULT: 1:3)
#' @param samps Targets data.frame - overrides subjects and frags
#'              (DEFAULT: NULL)
#' @param rightType "Correct" array type, tries to guess is NULL (DEFAULT: NULL)
#' @param force Parse IDAT files with different array size,
#'              see the man page ?read.metharray (DEFAULT: FALSE)
#' 
#' @return An RGChannelSet with pData $subject and $Basename filled 
#' 
#' @import  minfi
#' @import  sesame
#' 
#' @export 
#'
getRGChannelSet <- function(subjects = NULL,
                            frags = 1:3,
                            samps = NULL,
                            rightType = NULL,
                            force = FALSE) {

    if (is.null(samps)) {
        samps <- getSamps(subjects = subjects, frags = frags)
    }
    stopifnot(all(c("Basename", "subject") %in% names(samps)))

    rgSet <- readMethArrays(base = ".", targets = samps, verbose = TRUE,
                            rightType = rightType, force = force)
    sampleNames(rgSet) <- rgSet$subject

    return(rgSet)
}
