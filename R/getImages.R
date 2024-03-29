#' Get Paired Microscopic Images
#'
#' Return a list of paired DNA and fluorescent images
#' appropriate for \code{\link{parseImages}}.
#'
#' @param source A character vector identifying the source of images.
#'   The source can be a directory with subdirectories, each of which 
#'   contains images organized as specified by \code{which.images}. In 
#'   typical use, each subdirectory is the name of the well from a multi-
#'   well dish suchg as A1, C03, d003, etc. Alternatively, the source can be 
#'   a character string of a \code{.zip} file or character vector of image
#'   files organized as indicated by \code{which.images}. 
#' @param type A character string identifying the type of image files to parse
#'   ("tif", tiff", jpeg", jpg or png").
#' @param which.images An integer of length 2 or 3 or \code{NULL}. The first 
#'   number indicates the position of the DNA image. The second number 
#'   indicates the position of "target" image. The optional third number 
#'   specifies the total number of images for each field. If this is not 
#'   specified, the maxmimum of \code{which.images[1:2]} will be used 
#'   for this value. If \code{NULL}, no order is assumed and no consistency
#'   checks are performed. The default of \code{c(1, 2)} indicates a DNA image 
#'   followed by target image. A value of \code{c(2, 1)} indicates that 
#'   the target image is followed by the DNA image in pairs of images. A 
#'   value of \code{c(1, 2, 4)} indicates a DNA image, a target image, 
#'   and two additional images, which are ignored, in each set of four 
#'   images. 
#' @param pattern Optional grep pattern as character string used by
#'   \code{\link{list.files}} to select image files.
#' @param verbose If \code{TRUE}, print diagnostic messages as files are read.
#'
#' @details
#'
#' Images specified in \code{source} will be evaluated with the same logic in
#' \code{\link{checkImages}} to determine if the proper number of files are
#' present. The order required by \code{\link{parseImages}} is the DNA
#' image first and the fluorescent viral target second. The position of
#' the DNA image and the target image within each set of images can be changed
#' with the \code{which.images} argument.
#'
#' Images associated with each multiplicity of infection can be individual
#' files in a single directory where each directory named as the well such as
#' \code{A1}, \code{A2}, etc. and the files within are identified as
#' \code{A1/file001.tif}, \code{A1/file002.tif}, etc. The well identifier
#' can be in upper or lower case and can contain leading zeros such as
#' \code{c0003/file12.tif} as well as a numeric prefix such as \code{1A3}.
#'
#' Alternatively, each group of images associated with a given moi can be
#' a multi-layered tiff file where the sequence of images in the file is
#' specified by the argument \code{which.images}.
#'
#' If \code{'source'} is one or more zip files, files in the temporary directory 
#' (\code{\link{tempdir}}) will be deleted in order to receive the compressed
#' files. 
#'
#' @return
#'
#' A list of two named lists containing the nuclear and target images.
#' Each well or file will be represented as an element in the lists \code{nuc}
#' and \code{tgt}. Diagnostic messages are provided if \code{verbose}
#' is \code{TRUE}. 
#'
#' @examples
#' # Example with data organized by folder or well
#'   path.by.folder <- system.file("extdata", "by_folder", package = "virustiter")
#'   getImages(path.by.folder)
#'
#' @import EBImage
#'
#' @export
#'
getImages <- function(source, type = "tiff", which.images = c(1, 2, 2),
	pattern = NULL, verbose = FALSE)
{
# requires EBImage, ensure appropriate values for parameters
	if (!require(EBImage))
		stop("The 'EBImage' package must be installed with biocLite")

# verify source files or directories
	if (length(source) == 1 && !file.exists(source))
		stop("unable to find '", deparse(substitute(source)), "'")
	if (length(source) > 1 && !all(file.exists(source)))
		stop("not all files named in '", deparse(substitute(source)), "' exist")

# verify and adjust 'which.images' argument
	if (is.null(which.images))
		stop("Use 'readImage' directly to read a series of images")
	if (length(which.images) == 2)
			which.images <- c(which.images, max(which.images))
	if (length(which.images) != 3)
		stop("'which.images' must be an integer vector of length 2 or 3")
	if (which.images[3] != max(which.images))
		stop("the third value in 'which.images' must be the largest")

# collect image files, empty tempdir() for zip files
	if (all(file.info(source)$isdir)) # directory name(s)
		ff <- list.images(path = source, type = type, pattern = pattern)
	else if (all(grepl("zip$", source, ignore.case = TRUE))) { # zip file
		file.remove(list.files(tempdir(), full = TRUE, recursive = TRUE))
		unzip(source, exdir = tempdir())
		ff <- list.images(path = tempdir(), type = type, pattern = pattern)
	}
	else if (all(!file.info(source)$isdir)) # file name(s)
		ff <- source
	else
		stop("unable to use files/source in ", deparse(substitute(source)))
	if (verbose) {
		message("Found ", length(ff), " image file", ifelse(length(ff) == 1, "", "s"))
		flush.console()
	}

# extract fields to determine if images are organized by well or stack
	spl <- strsplit(ff, "/")
	field1 <- sapply(spl, tail, 1) # last field, file name
	field2 <- sapply(spl, function(x) head(tail(x, 2), 1)) # potential well name
	wellpat <- "[[:alpha:]][[:digit:]]+$" # pattern for 'well' at end of string
	sel <- grepl(wellpat, field2)

# assign value to imageType as "byWell" or "byFile" and complete message
	if (all(sel)) { # extract well and numeric optional prefix
		imageType <- "byWell"
		plate <- well.info(field2)$plate
		well <- well.info(field2)$well
		filename <- NULL
	}
	else if (!any(sel)) {
		imageType <- "byFile"
		well <- NULL
		filename <- field1
	}
	else
		stop("unable to use mixture of image files in ",
			deparse(substitute(source)), '"')

# split image paths into related groups (by field2 or by file)
	if (imageType == "byWell")
		ffsplit <- split(ff, field2)
	else if (imageType == "byFile")
		ffsplit <- split(ff, filename)
	else
		stop("CAN'T HAPPEN! Unexpected value for 'imageType'")

# read all images as a list and coerce to grayscale with a warning
	if (verbose) {
		message("Reading images grouped by ", ifelse(imageType == "byWell", "well", "file"))
		flush.console()
	}
	img <- lapply(ffsplit, function(f) suppressWarnings(readImage(f)))
	if (any(sapply(img, colorMode) != 0)) {
		warning("images have been converted to grayscale by uniform RGB averaging",
			call. = FALSE)
		img <- lapply(img, channel, "gray")
	}

# perform check for expected numbers of images
	n_dna <- which.images[1]
	n_tgt <- which.images[2]
	n_field <- which.images[3]

# are total images in each group sensible?
	n <- sapply(img, numberOfFrames)
	bad <- n < n_dna | n < n_tgt | n < n_field
	if (all(bad))
		stop("All ", sum(bad), " images have fewer frames than specified", 
		" in 'which.images'")
	if (any(bad))
		stop(sum(bad), " of ", length(bad),
			"had fewer frames than specified in 'which.images'")

# are total images in each group a multiple of field size?
	bad <- which(n %% n_field != 0)
	if (length(bad > 10))
		warning("many images in ", deparse(substitute(source)),
			" are not multiples of ", n_field, call. = FALSE)
	else if (length(bad > 1))
		warning("the number of images in:\n", paste(names(img)[bad], collapse = ", "),
			" are not multiples of ", n_field, call. = FALSE)
	else if (length(bad == 1))
		warning("the number of images in ", deparse(substitute(source)),
			" is not a multiple of ", n_field, call. = FALSE)

	# extract dna and target images
		idx <- lapply(n, function(N) seq(n_dna, N, n_field))
		dnaImages <- Map(function(x, i) x[,,i], img, idx)

		idx <- lapply(n, function(N) seq(n_tgt, N, n_field))
		tgtImages <- Map(function(x, i) x[,,i], img, idx)

# count and report on the number of groups to display
	nspl <- length(ffsplit)
	msg <- sprintf("Found %d group%s of image %s", nspl,
		ifelse(nspl == 1, "", "s"),
		switch(n_field, "singles", "pairs", "triplets", "quads", "sets", "sets"))
	if (verbose) {
		message(msg)
		flush.console()
	}

# report on findings 
	if (verbose) {
		id <- names(ffsplit)
		count <- sapply(dnaImages, function(v) dim(v)[3])
		count[is.na(count)] <- 1
		message(sprintf("%3d: %d pair%s in %s\n", seq_along(id), count,
			ifelse(count == 1, "", "s"), id), appendLF = FALSE)
		flush.console()
	}

# return image pairs
	return(list(nuc = dnaImages, tgt = tgtImages))
}
