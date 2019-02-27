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
#' @param pattern Optional grep pattern as character string used by \code{\link{list.files}}
#'   to select image files.
#' @param method Character string specifying the method of displaying images.
#'   Default of \code{"none"} simply summarizes the images. A value of \code{"raster"}
#'   uses R raster graphics and \code{"browser"} attempts to use a browser.
#'   (Unfortunately this seems to be failing with EBImage version 4.22.0.)
#' @param ask Logical value to use \code{par(ask = TRUE)} if \code{method = "raster"}.
#' @param separate Logical value to normalize each frame separately if
#'   \code{TRUE}. Aplies only if \code{display} is not \code{"none"}. 
#' @param nx Integer value passed to \code{display{}} that specifies the number of
#'   images in a column if \code{display} is not \code{"none"}.
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
#' If \code{'source'} is a zip file, files in the temporary directory 
#' (\code{\link{tempdir}}) will be deleted in order to receive the compressed
#' files. 
#'
#' @return
#'
#' An \code{invisible} character vector of image files, diagnostic messages
#' are printed on the console with the option to display normalized image pairs.
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
	pattern = NULL, method = c("none", "raster", "browser"), ask = TRUE,
	separate = FALSE, nx = 2)
{
# requires EBImage, ensure appropriate values for parameters
	if (!require(EBImage))
		stop("The 'EBImage' package must be installed with biocLite")
	method <- match.arg(method)

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

# assign 'method'
	if (method == "browser") {
		message("Currently the display with browser seems to fail in EBImage")
		message("The display method has been changed to 'raster'")
		flush.console()
		method <- "raster"
		opar <- par(ask = TRUE)
	}

# collect image files, empty tempdir() for zip files
	if (all(file.info(source)$isdir))
		ff <- list.images(path = source, type = type, pattern = pattern)
	else if (all(grepl("zip$", source, ignore.case = TRUE))) {
		file.remove(list.files(tempdir(), full = TRUE, recursive = TRUE))
		unzip(source, exdir = tempdir())
		ff <- list.images(path = tempdir(), type = type, pattern = pattern)
	}
	else if (all(!file.info(source)$isdir))
		ff <- source # must be character vector of file names
	else
		stop("unable to use files/source in ", deparse(substitute(source)))

# provide status messages
	nff <- length(ff)
	txt <- paste0("Found ", nff, " image file", ifelse(nff == 1, " ", "s "))
	message(txt, appendLF = FALSE)
	flush.console()

# extract fields to determine if images are organized by well or stack
	spl <- strsplit(ff, "/")
	field1 <- sapply(spl, tail, 1)
	field2 <- sapply(spl, function(x) head(tail(x, 2), 1))
	sel <- grepl("^[abcdefghijklmnop][[:digit:]]+$", field2, ignore.case = TRUE)

# assign value to imageType as "byWell" or "byFile" and complete message
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
		stop("unable to use mixture of image files in ",
			deparse(substitute(source)), '"')
	message("grouped by ", ifelse(imageType == "byWell", "well", "file"))

# split image paths into related groups (by well or by file)
	if (imageType == "byWell")
		ffsplit <- split(ff, well)
	else if (imageType == "byFile")
		ffsplit <- split(ff, filename)
	else
		stop("can't happen! Unexpected value for 'imageType'")

# read all images as a list and coerce to grayscale with a warning
	message("Reading images...", appendLF = FALSE)
	flush.console()
	img <- lapply(ffsplit, function(f) suppressWarnings(readImage(f)))
	message("done")
	flush.console()
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

# extract dna images and adjust to 3 dimensions to use "["
	dnaImages <- lapply(img, function(x, first = n_dna, by = n_field) {
		dm <- dim(x)
		if (length(dm) == 2)
			dim(x) <- c(dm, 1)
		N <- dim(x)[3]
		x[,,seq(first, N, by), drop = TRUE]})	

# extract tgt images and adjust to 3 dimensions to use "["
	tgtImages <- lapply(img, function(x, first = n_mfi, by = n_field) {
		dm <- dim(x)
		if (length(dm) == 2)
			dim(x) <- c(dm, 1)
		N <- dim(x)[3]
		x[,,seq(first, N, by), drop = TRUE]})

# count and report on the number of groups to display
	nspl <- length(ffsplit)
	message("Found ", nspl, " group", ifelse(nspl == 1, "", "s"), " of image pairs")
	flush.console()

# report on and process each group of files in turn
	for (IDX in seq_along(ffsplit)) {
		nimg <- dim(dnaImages[[IDX]])[3]
		if(is.na(nimg)) nimg <- 1
		message(sprintf("%3d: %d in ", IDX, nimg), names(ffsplit)[IDX])
		if (method != "none") {
			myDna <- normalize(dnaImages[[IDX]], separate = separate)
			myTgt <- normalize(tgtImages[[IDX]], separate = separate)
			n <- dim(myDna)[3]
			if (is.null(n) || is.na(n)) n <- 1
			i <- c(rbind(seq_len(n), n + seq_len(n)))
			img <- combine(myDna, myTgt)[,,i]
			opar <- par(ask = ask)
			display(tile(img, nx = nx), all = TRUE, method = method,
				title = names(ffsplit)[IDX])
			par(opar)
		}
	}
	invisible(ff) # return list of file names
}

