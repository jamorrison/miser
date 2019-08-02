#' download IDAT file from a provided URL
#'
#' Checks to see if the requested IDAT file exists in the current working directory.
#' Downloads the IDAT if it does not exist. Returns basename of IDAT file.
#'
#' @param IDAT URL for IDAT file to attempt downloading
#'
#' @return     basename of IDAT file (see R basename function for more info)
#'
#' @export
#'
getIDAT <- function(IDAT) {
  fname <- basename(IDAT)
  if (fname %in% list.files()) {
    message("Found ", fname, " locally; skipping download.") 
  } else { 
    message("Downloading ", fname, " from ", IDAT, "...", appendLF=FALSE)
    download.file(IDAT, fname) 
    message("done.")
  }
  return(fname)
}
