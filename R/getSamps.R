#' Helper function to read in subjects (samples) and Basename columns
#' 
#' @param subjects The names of each subject, will autodetect if NULL
#' (DEFAULT: NULL)
#' @param frags Which elements to extract for the array Basename, ignored if
#' `check.frags` = TRUE (DEFAULT: 1:3)
#' @param check.frags Check if frags covers all elements up to "idat.gz"
#' (DEFAULT: TRUE)
#' 
#' @return A data.frame with Basename and subject columns
#' 
#' @export 
getSamps <- function(subjects=NULL, frags=1:3, check.frags=TRUE) { 
    # short circuit all of this if the subjects are a list of IDATs:
    if (all(grepl("idat", ignore.case=TRUE, subjects))) { 
        # list of IDATS
        if (check.frags) { # find last element before "idat.gz"
            subjects <- sub("_Grn.idat.gz", "", subjects)
            subjects <- sub("_Red.idat.gz", "", subjects)
            samps <- data.frame(Basename = unique(subjects),
                                subject = unique(elts(subjects)),
                                stringsAsFactors = FALSE)
            return(samps)
        } else {
            samps <- data.frame(Basename=unique(elts(subjects, elt=frags)),
                                subject=unique(elts(subjects)),
                                stringsAsFactors = FALSE)
            return(samps)
        }
    } else { 
        basenames <- unique(elts(list.files(patt="*idat*"), elt=frags))
        samps <- data.frame(Basename=basenames, stringsAsFactors = FALSE)
        if (is.null(subjects)) {
            samps$subject <- elts(samps$Basename)
        } else { 
            stopifnot(identical(names(subjects), samps$Basename))
        }
        return(samps)
    }
}
