#' Read Illumina methylation arrays
#'
#' Code based on minfi::read.metharray and minfi::read.metharray.exp. Additional
#' functionality was needed beyond the minfi function as it is possible to find
#' datasets which have been incorrectly labeled as 450k or EPIC arrays.
#'
#' @param base Base directory of IDAT files
#' @param targets data.frame of target Red/Green IDAT files
#' @param extended Return a RGChannelSetExtended
#' @param verbose Output extra information
#' @param force Parse IDAT files with different array size, see the man
#'              page ?read.metharray for more information
#' @param dropWrongType Drop "wrong" array types (DEFAULT: TRUE)
#' @param rightType "Correct" array type, tries to guess if NULL
#'                  (DEFAULT: NULL)
#'
#' @return A RGChannelSet
#'
#' @importFrom BiocParallel bplapply
#' @importFrom illuminaio readIDAT
#' @import minfi
#'
#' @seealso minfi::read.metharray() and minfi::read.metharray.exp()
#'
#' @examples
#'
#' @export
#'
readMethArrays <- function(base,
                           targets,
                           extended = FALSE,
                           verbose = FALSE,
                           force = FALSE,
                           dropWrongType = TRUE,
                           rightType = NULL) {

    # Setup up Bioconductor parallelization parameters
    # TODO: When Tim invariably asks for this to be HDF5-backed, this code will
    #       will need to be tested out to see if it is behaving properly. Until
    #       then, it will be left commented out and is unneeded.
    # TODO: Make BPPARAM an argument the user can set in the function call
    #BACKEND <- getRealizationBackend()
    BPREDO <- list()
    BPPARAM <- SerialParam()

    # Check the provided files exist
    G.files <- .checkFilesExist(base = base, targets = targets, color = "green")
    R.files <- .checkFilesExist(base = base, targets = targets, color = "red")

    # Load "Quants" from IDAT file, include "SD", and "NBeads" if
    # extended = TRUE
    stime <- system.time({
        G.Quants <- bplapply(G.files,
                             .readIDATQuants,
                             extended = extended,
                             verbose = verbose,
                             BPREDO = BPREDO,
                             BPPARAM = BPPARAM)
        R.Quants <- bplapply(R.files,
                             .readIDATQuants,
                             extended = extended,
                             verbose = verbose,
                             BPREDO = BPREDO,
                             BPPARAM = BPPARAM)
    })[3]

    if (verbose) {
        message(sprintf("[readMethArrays] Read IDAT files in %.1f seconds",
                        stime))
        message("[readMethArrays] Creating data matrices...", appendLF = FALSE)
    }

    # Check to make sure files are of the same length and array type
    ptime1 <- proc.time()
    allNProbes <- vapply(G.Quants, nrow, integer(1L))
    arrayTypes <- cbind(do.call(rbind, lapply(allNProbes, .guessArrayType)),
                        size = allNProbes)
    sameLength <- (length(unique(arrayTypes[, "size"])) == 1)
    sameArray <- (length(unique(arrayTypes[, "array"])) == 1)

    if (!sameLength && !sameArray && !dropWrongType) {
        message("[readMethArrays] Trying to parse IDAT files from different ",
                "arrays...")
        message("    Guessed array sizes and types:")
        print(arrayTypes[, c("array", "size")])
        stop("[readMethArrays] Trying to parse IDAT files with different ",
             "sizes and types")
    } else if (!sameLength && !sameArray && dropWrongType) {
        if (is.null(rightType)) {
            rightType <- .pickMostCommonArrayType(arrayTypes)
        } else {
            rightType <- .getArrayType(rightType)
        }

        message("[readMethArrays] Trying to parse IDAT files from different ",
                "arrays...")
        message("[readMethArrays] 'dropWrongType' == TRUE, so arrays not of ",
                "type ", rightType, " will be dropped from data")

        G.Quants <- subset(G.Quants, arrayTypes == rightType)
        R.Quants <- subset(R.Quants, arrayTypes == rightType)
    }
    if (!sameLength && sameArray && !force) {
        stop("[readMethArrays] Trying to parse IDAT files with different ",
             "array size. However, they do appear to all have the same array ",
             "type. You can force this by setting 'force = TRUE'. See man ",
             "page ?readMethArrays for more information.")
    }

    commonAddresses <- as.character(Reduce("intersect", lapply(G.Quants,
                                                               rownames)))

    # NOTE: colnames must be manually set as it is not safe to assume they will
    #       be correctly parsed
    GreenMean <- do.call(cbind,
                         lapply(G.Quants,
                                function(xx) xx[commonAddresses, "Mean"]))
    colnames(GreenMean) <- names(G.Quants)
    RedMean <- do.call(cbind,
                       lapply(R.Quants,
                              function(xx) xx[commonAddresses, "Mean"]))
    colnames(RedMean) <- names(R.Quants)

    if (extended) {
        GreenSD <- do.call(cbind,
                           lapply(G.Quants,
                                  function(xx) xx[commonAddresses, "SD"]))
        colnames(GreenSD) <- names(G.Quants)
        RedSD <- do.call(cbind,
                         lapply(R.Quants,
                                function(xx) xx[commonAddresses, "SD"]))
        colnames(RedSD) <- names(R.Quants)
        NBeads <- do.call(cbind,
                          lapply(G.Quants,
                                 function(xx) xx[commonAddresses, "NBeads"]))
        colnames(NBeads) <- names(G.Quants)
    }

    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        message(sprintf("done in %.1f seconds", stime))
        message("[readMethArrays] Instantiating final object...",
                appendLF = FALSE)
    }

    ptime1 <- proc.time()
    if (extended) {
        out <- RGChannelSetExtended(Red = RedMean,
                                    Green = GreenMean,
                                    RedSD = RedSD,
                                    GreenSD = GreenSD,
                                    NBeads = NBeads)
    } else {
        out <- RGChannelSet(Red = RedMean, Green = GreenMean)
    }
    rownames(out) <- commonAddresses

    out@annotation <- c(array = arrayTypes[1, 1], annotation = arrayTypes[1, 2])

    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if (verbose) {
        message(sprintf("done in %.1f seconds", stime))
    }

    pD <- subset(targets, Basename %in% colnames(out))
    #pD <- targets
    pD$filenames <- file.path(base, colnames(out))
    rownames(pD) <- colnames(out)
    colData(out) <- as(pD, "DataFrame")

    return(out)
}

### Check if files exist, stops if 1 or more files don't exist
###
### PARAMETERS
### base    ---> directory path to targets
### targets ---> data.frame with 'Basename' of individual IDATs
### color   ---> which color (red or green) of IDATs to check
###
### RETURNS
### character vector of files existing files
.checkFilesExist <- function(base,
                             targets,
                             color = c("red", "green")) {
    # Check that the target basenames exist in the targets data.frame
    if (!"Basename" %in% names(targets)) {
        stop("Need 'Basename' among the column names of 'targets'")
    }

    # Set which color (green or red) IDATs to check
    color <- match.arg(color)
    if (color == "green") {
        suffix <- "_Grn.idat"
    } else if (color == "red") {
        suffix <- "_Red.idat"
    }

    files <- file.path(base, basename(targets$Basename))

    files <- sub("_Grn\\.idat.*", "", files)
    files <- sub("_Red\\.idat.*", "", files)
    stopifnot(!anyDuplicated(files))

    out.files <- paste0(files, suffix)
    names(out.files) <- basename(files)
    these.dont.exist <- !file.exists(out.files)
    if (any(these.dont.exist)) {
        out.files[these.dont.exist] <- paste0(out.files[these.dont.exist],
                                              ".gz")
    }
    if (!all(file.exists(out.files))) {
        not.found <- sub("\\.gz", "", out.files[!file.exists(out.files)])
        stop("The following files do not exist: ",
             paste(not.found, collapse = ", "))
    }

    return(out.files)
}

### Read "Quants" entry from IDAT file
###
### PARAMETERS
### files    ---> IDAT files to read
### extended ---> Keep "SD" and "NBeads" for extended RGChannelSet
### verbose  ---> Print out messages?
###
### RETURNS
### matrix of the average intensity ("Mean", always returned)
###               number of beads ("NBeads", returned if extended = TRUE)
###               measure of variability ("SD", returned if extended = TRUE)
.readIDATQuants <- function(files,
                            extended = FALSE,
                            verbose = FALSE) {
    # TODO: When BACKEND is added in, make sure to include BACKEND = NULL inside
    #       this function definition (right after verbose). Will also need to
    #       add this to the function description above.
    if (verbose) {
        message("[readMethArrays] Reading ", basename(files))
    }

    Quants <- illuminaio::readIDAT(files)[["Quants"]]
    if (!extended) {
        Quants <- Quants[, "Mean", drop = FALSE]
    }
    # TODO: Uncomment when BACKEND has been added in
    #if (!is.null(BACKEND)) {
    #    Quants <- realize(Quants, BACKEND = BACKEND)
    #}

    return(Quants)
}

### Try to guess the type of array based on the number of probes
###
### PARAMETERS
### nProbes ---> Number of probes in file
###
### RETURNS
### named character vector of predicted array and annotation
.guessArrayType <- function(nProbes) {
    .default.27k.annotation <- "ilmn12.hg19"
    .default.450k.annotation <- "ilmn12.hg19"
    .default.epic.annotation <- "ilm10b4.hg19"

    if (nProbes >= 41000 & nProbes <= 41100) {
        # Unknown methylation chip
        arrayAnnotation <- c(array = "HorvathMammalMethylChip40",
                             annotation = "test.unknown")
    } else if (nProbes >= 54000 & nProbes <= 56000) {
        # Illumina 27k methylation array
        arrayAnnotation <- c(array = "IlluminaHumanMethylation27k",
                             annotation = .default.27k.annotation)
    } else if (nProbes >= 622000 & nProbes <= 623000) {
        # Illumina 450k methylation array
        arrayAnnotation <- c(array = "IlluminaHumanMethylation450k",
                             annotation = .default.450k.annotation)
    } else if (nProbes >= 1032000 & nProbes <= 1033000) {
        # Illumina EPIC (850k) methylation array (old scan type)
        arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC",
                             annotation = .default.epic.annotation)
    } else if (nProbes >= 1050000 & nProbes <= 1053000) {
        # Illumina EPIC (850k) methylation array (current scan type)
        arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC",
                             annotation = .default.epic.annotation)
    } else {
        arrayAnnotation <- c(array = "Unknown",
                             annotation = "Unknown")
    }

    return(arrayAnnotation)    
}

### Pick array type with the greatest number of occurrences
###
### PARAMETERS
### arrays ---> matrix containing the array types
###
### RETURNS
### name of the array with the greatest number of occurrences
.pickMostCommonArrayType <- function(arrays) {
    # TODO: Only select mostCommon if it has a large percentage of the number of
    #       arrays (say >75%)
    counts <- as.data.frame(table(arrays[,"array"]))
    mostCommon <- as.character(counts$Var1[which.max(counts$Freq)])

    return(mostCommon)
}

### Retrieve the name of the array type based on input value
###
### PARAMETERS
### arrayName ---> either a character or numeric value of the array type
###                     valid character values: "27k", "450k", "EPIC"
###                     valid numeric values: 8490, 13534, 21145
###                 numeric values correspond to the GEO platform numbers
###
### RETURNS
### full name of the input array name
.getArrayType <- function(arrayName) {
    if (is.character(arrayName)) {
        if (arrayName == "27k") {
            type <- "IlluminaHumanMethylation27k"
        } else if (arrayName == "450k") {
            type <- "IlluminaHumanMethylation450k"
        } else if (arrayName == "EPIC") {
            type <- "IlluminaHumanMethylationEPIC"
        } else {
            stop("Unknown array type (", arrayName, ") provided.")
        }
    } else if (is.numeric(arrayName)) {
        if (arrayName == 8490) {
            type <- "IlluminaHumanMethylation27k"
        } else if (arrayName == 13534) {
            type <- "IlluminaHumanMethylation450k"
        } else if (arrayName == 21145) {
            type <- "IlluminaHumanMethylationEPIC"
        } else {
            stop("Unknown array type (", arrayName, ") provided.")
        }
    } else {
        stop("Unknown class of array type (", class(arrayName), ") provided.")
    }

    return(type)
}
