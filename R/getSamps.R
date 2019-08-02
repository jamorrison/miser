#' helper function to read in samples and Basename columns
#' 
#' @param   subjects     the names of each subject, or NULL to autodetect (NULL)
#' @param   frags        which elements to extract for the array Basename (1:3)
#' @param   check.frags  check if frags covers all elements up to "idat.gz" (FALSE)
#' 
#' @return            a data.frame
#' 
#' @export 
getSamps <- function(subjects=NULL, frags=1:3, check.frags=FALSE) { 
  # short circuit all of this if the subjects are a list of IDATs:
  if (all(grepl("idat", ignore.case=TRUE, subjects))) { 
    # list of IDATS
    if (check.frags) { # find last element before "idat.gz"
      idx <- sapply(strsplit(subjects,"_",fixed=TRUE), function(x) grep("idat.gz$",x))
      stopifnot(all(abs(idx - mean(idx)) < 10^5)) # check all entries have the same number of elements
      frags <- 1:(idx[1]-1) #idx is the index with "idat.gz" in it, want to stop at element before this
    }
    samps <- data.frame(Basename=unique(elts(subjects, elt=frags)),
                        subject=unique(elts(subjects)))
  } else { 
    basenames <- unique(elts(list.files(patt="*idat*"), elt=frags))
    samps <- data.frame(Basename=basenames)
    if (is.null(subjects)) {
      samps$subject <- elts(samps$Basename)
    } else { 
      stopifnot(identical(names(subjects), samps$Basename))
    }
  }
  return(samps)
}
