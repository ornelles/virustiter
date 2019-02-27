#' Get Paired Microscopic Images
#'
#' Return a list of paired DNA and fluorescent images
#' appropriate for \code{\link{parseImages}}.
#'
#' @param source A character vector identifying a directory or directories
#'   with multilayer tiff files \emph{or} subdirectories identified
#'   by well with separate, paired images \emph{or} a character vector
#'   of image files \emph{or} a \code{.zip} file with the above.
#' @param type A character string identifying the type of image files to parse
#'   ("tif", "tiff", "jpeg", "jpg" or "png".)
#' @param which.images An integer of length 2 or 3. The first two numbers indicate
#'   the relative position of the DNA image and the target in each field. The optional
#'   third number specifies the total number of images for each field. A value of
#'   \code{c(1, 2)} indicates DNA first and target second. A value of \code{c(2, 1)}
#'   indicates that the order is target first and DNA image second. A value of
#'   \code{c(1, 2, 3)} indicates a DNA image, a target image, and a third (ignored)
#'   image such as a phase contrast image or second fluorescent color in each set.
#' @param pattern Optional character string to serve as a \code{grep} pattern
#'   for \code{\link{list.files}} to select image files.
#' @param verbose If \code{TRUE}, print diagnostic messages as files are read.
#'
#' @details
#'
#' Images specified in \code{source} will be checked for the proper
#' organization of image files. The order required by \code{\link{parseImages}}
#' is the DNA image first and fluorescent viral target second. This
#' order can be changed with the \code{which.images} argument.
#'
#' Images associated with each multiplicity of infection can be individual
#' files in a single directory where each directory named as the well such as
#' \code{A1}, \code{A2}, etc. and the files within are identified as
#' \code{A1/file001.tif}, \code{A1/file002.tif}, etc. The well identifier
#' can be in upper or lower case and can contain leading zeros such as
#' \code{c0003/file12.tif}.
#'
#' Alternatively, each group of images associated with a given moi can be
#' a multi-layered tiff file where the sequence of images in the file is
#' specified by the argument \code{which.images}.
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
	if (verbose)
		message("Found ", length(ff), " image file", ifelse(length(ff) == 1, "", "s"))

# extract fields to determine if images are organized by well or stack
	spl <- strsplit(ff, "/")
	field1 <- sapply(spl, tail, 1)
	field2 <- sapply(spl, function(x) head(tail(x, 2), 1))
	sel <- grepl("^[abcdefghijklmnop][[:digit:]]+$", field2, ignore.case = TRUE)

# assign variables to direct processing
	if (all(sel)) {
		imageType <- "byWell"
		well <- field2
		filename <- NULL
	}
	else if (!any(sel)) {
		imageType <- "byFile"
		well <- NULL
		filename <- field1
	}
	else
		stop("unable to use mixture of image files in ", source, '"')

# split image paths into related groups (by well or by file)
	if (imageType == "byWell")
		ffsplit <- split(ff, well)
	else if (imageType == "byFile")
		ffsplit <- split(ff, filename)
	else
		stop("internal error: unexpected value for 'imageType'")

# read all images as a list and coerce to grayscale with a warning
	if (verbose)
		message("Reading images grouped by ", ifelse(imageType == "byWell", "well", "file"))
	img <- lapply(ffsplit, function(f) suppressWarnings(readImage(f)))
	if (any(sapply(img, colorMode) != 0)) {
		warning("images have been converted to grayscale by uniform RGB averaging")
		img <- lapply(img, channel, "gray")
	}

# check for expected numbers of images
	n_dna <- which.images[1]
	n_mfi <- which.images[2]
	n_field <- which.images[3]

# are total images in each group sensible?
	n <- sapply(img, function(x) dim(x)[3])
	bad <- which(n %% n_field != 0)
	if (length(bad > 1))
		stop("\nThe number of images in ", paste(names(img)[bad], collapse = ", "),
			" are not multiples of ", n_field)
	else if (length(bad == 1))
		stop("\nThe number of images in ", deparse(substitute(source)), 
				" is not a multiple of ", n_field)

# extract dna images
	dnaImages <- lapply(img, function(x, first = n_dna, by = n_field) {
		dm <- dim(x)
		if (length(dm) == 2)
			dim(x) <- c(dm, 1)
		N <- dim(x)[3]
		x[,,seq(first, N, by), drop = TRUE]})

# extract mfi images and adjust to 3 dimensions
	mfiImages <- lapply(img, function(x, first = n_mfi, by = n_field) {
		dm <- dim(x)
		if (length(dm) == 2)
			dim(x) <- c(dm, 1)
		N <- dim(x)[3]
		x[,,seq(first, N, by), drop = TRUE]})

# count and report on the number of groups to display
	nff <- length(ffsplit)
	if (verbose)
		message("Found ", nff, " group", ifelse(nff == 1, "", "s"), " of images")

# report on findings 
	if (verbose) {
		id <- names(ffsplit)
		count <- sapply(dnaImages, function(v) dim(v)[3])
		count[is.na(count)] <- 1
		message(sprintf("%3d: %d pair%s in %s\n", seq_along(id), count,
			ifelse(count == 1, "", "s"), id), appendLF = FALSE)
	}

# return image pairs
	return(list(nuc = dnaImages, tgt = mfiImages))
}
