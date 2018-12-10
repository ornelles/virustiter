#' Check Paired Microscopic Images
#'
#' Check validity of paired DNA and a fluorescent images
#' appropriate for \code{parseImages()}.
#'
#' @param source A character vector identifying a directory or directories
#'   with multilayer tiff files \emph{or} subdirectories identified
#'   by well with separate, paired images per well \emph{or} a character vector
#'   of image files \emph{or} a \code{.zip} file with the above.
#' @param type A character string identifying the type of image files to parse
#'   ("tif", "tiff", "jpeg", "jpg" or "png".)
#' @param which.images An integer of length 2 or 3. The first two numbers indicate
#'   the relative position of the DNA image and the target image in each field. The optional
#'   third number specifies the total number of images for each field. A value of
#'   \code{c(1, 2)} indicates DNA first and target second. A value of \code{c(2, 1)}
#'   indicates that the order is target first and DNA image second. A value of
#'   \code{c(1, 2, 3)} indicates a DNA image, a target image, and a third (ignored)
#'   image such as a phase contrast image or second fluorescent color in each set.
#' @param pattern Optional grep pattern as character string used by \code{list.files()}
#'   to select image files.
#' @param method Character string specifying the method of displaying images.
#'   Default of \code{"none"} simply summarizes the images. A value of \code{"raster"}
#'   uses R raster graphics and \code{"browser"} attempts to use a browser.
#'   (Unfortunately this seems to be failing with EBImage version 4.22.0.)
#' @param ask Logical value to use \code{par(ask = TRUE)} if \code{method = "raster"}.
#'
#' @details
#'
#' Images specified in \code{source} will be evaluated with the same logic in
#' \code{\link{getImages}} to determine if the proper number of files are
#' present and report on the number and form of the image files. 
#' The default value of the \code{which.images} argument treats the first of each
#' pair a DNA image and the second as a fluorescent 
#' image of the viral target.
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
#' This function is called for the side effect of examining the sequence
#' of image files and optionally displaying the normalized images.
#' Diagnostic messages are printed on the console with \code{message()}.
#' Use the \code{grep} string in \code{pattern} to limit the selected images.
#' A character vector of the selected image files will invisibly returned. 
#'
#' @examples
#' # Example with data organized by folder or well
#'   path.by.folder <- system.file("extdata", "by_folder", package = "virustiter")
#'   checkImages(path.by.folder)
#'
#' @import EBImage
#'
#' @export
#'
checkImages <- function(source, type = "tiff", which.images = c(1, 2, 2),
	pattern = NULL, method = c("none", "raster", "browser"),
	ask = ifelse(method == "raster", TRUE, NULL))
{
# requires EBImage, ensure appropriate values for parameters
	if (!require(EBImage))
		stop("The 'EBImage' package must be installed with biocLite")

# are all files or directories found in 'source' argument legitimate?
	if (!all(file.exists(source)))
		stop("not all files named in ", deparse(substitute(source)), " exist")

# collect image files
	if (all(file.info(source)$isdir))
		ff <- list.images(path = source, type = type, pattern = pattern)
	else if (all(grepl("zip$", source, ignore.case = TRUE))) {
		unzip(source, exdir = tempdir())
		ff <- list.images(path = tempdir(), type = type, pattern = pattern)
	}
	else if (all(!file.info(source)$isdir))
		ff <- source
	else
		stop("unable to use files/source in ", deparse(substitute(source)))
	message("Found ", length(ff), " image files"); flush.console()

# check on arguments
	if (length(which.images) == 2)
			which.images <- c(which.images, max(which.images))
	if (length(which.images) != 3)
		stop("'which.images' must be an integer vector of length 2 or 3")
	if (which.images[3] != max(which.images))
		stop("the third value in 'which.images' must be the largest")

# method
	method <- match.arg(method)
	if (method == "browser") {
		message("Currently the display with browser seems to fail in EBImage")
		message("The value has been changed to 'raster'.")
		message("Use 'par(ask = TRUE)' to view images one-by-one")
		flush.console()
		method <- "raster"
	}

# extract fields to determine if images are organized by well or stack
	spl <- strsplit(ff, "/")
	field1 <- sapply(spl, tail, 1)
	field2 <- sapply(spl, function(x) head(tail(x, 2), 1))
	sel <- grepl("[[:alpha:]][[:digit:]]+$", field2) # test for well pattern

# assign variables to direct processing
	if (all(sel)) {
		imageType <- "byWell"
		well <- field2
		filename <- NULL
	}
	else if (!any(sel)) {
		imageType <- "byStack"
		well <- NULL
		filename <- field1
	}
	else
		stop("unable to use mixture of image files in ", source, '"')

# split image paths into related groups (by well or by file)
	if (imageType == "byWell")
		ffsplit <- split(ff, well)
	else if (imageType == "byStack")
		ffsplit <- split(ff, filename)
	else
		stop("unexpected value for 'imageType'")

# read all images as a list and coerce to grayscale with a warning
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
		stop("\nThe number of images in ",
			" is not a multiple of ", n_field)

# extract dna images and adjust to 3 dimensions
	dnaImages <- lapply(img, function(x, first = n_dna, by = n_field) {
		dm <- dim(x)
		if (length(dm) == 2)
			dim(x) <- c(dm, 1)
		N <- dim(x)[3]
		x[,,seq(first, N, by)]})

# extract mfi images and adjust to 3 dimensions
	mfiImages <- lapply(img, function(x, first = n_mfi, by = n_field) {
		dm <- dim(x)
		if (length(dm) == 2)
			dim(x) <- c(dm, 1)
		N <- dim(x)[3]
		x[,,seq(first, N, by)]})

# count and report on the number of groups to display
	nff <- length(ffsplit)
	message("Found ", nff, " groups of images"); flush.console()

# process each group of files in turn
	for (IDX in seq_along(ffsplit)) {
		message(sprintf("%4d: ", IDX), "  ", names(ffsplit)[IDX])
		myDna <- normalize(dnaImages[[IDX]])
		myMfi <- normalize(mfiImages[[IDX]])
		n <- dim(myDna)[3]
		if (is.null(n) || is.na(n))
			n <- 1
		i <- c(rbind(seq_len(n), n + seq_len(n)))
		img <- combine(myDna, myMfi)[,,i]
		if (method != "none")
			opar <- par(ask = ask)
			display(tile(img, nx = 2), all = TRUE, method = method,
				title = names(ffsplit)[IDX])
			par(opar)
	}
	message("Done")
	invisible(ff)
}
