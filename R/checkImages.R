#' Check Paired Microscopic Images
#'
#' Check validity of paired DNA and a fluorescent images
#' appropriate for \code{parseImages()}.
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
#' @param method Character string specifying the method of displaying images.
#'   Default of \code{"none"} simply summarizes the images. A value of \code{"raster"}
#'   uses R raster graphics and \code{"browser"} attempts to use a browser.
#'   (Unfortunately this seems to be failing with EBImage version 4.22.0.)
#' @param ask Logical value to use \code{par(ask = TRUE)} if \code{method = "raster"}.
#' @param separate Logical value to normalize each frame separately if
#'   \code{TRUE}. Aplies only if \code{method} is \code{"raster"} or \code{"browser"}. 
#' @param nx Integer value passed to the \code{display} function that specifies
#'   the number of images in a column if \code{method} is not \code{"none"}.
#'
#' @details
#'
#' If \code{which.images} is not \code{NULL}, the images specified in
#' \code{source} will be evaluated with the same logic in
#' \code{\link{getImages}} to determine if the proper number of files are
#' present and report on the number and form of the image files. If
#' \code{which.images} is \code{NULL}, no check will be performed. 
#'
#' Images associated with each multiplicity of infection can be individual
#' files in a single directory where each directory is named for the well
#' such as \code{A1}, \code{A2}, etc. and the files within are identified as
#' \code{A1/file001.tif}, \code{A1/file002.tif}, etc. The well identifier
#' can be in upper or lower case and can contain leading zeros such as
#' \code{c0003/file12.tif}/ The well identifier also can contain a leading
#' numeric prefix such as \code{1A2} or \code{02H012}.
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

# provide immediate warnings
	owarn <- options(warn = 1)
	on.exit(options(owarn))
	
# verify source files or directories
	if (length(source) == 1 && !file.exists(source))
		stop("unable to find '", deparse(substitute(source)), "'")
	if (length(source) > 1 && !all(file.exists(source)))
		warning("not all files named in '", deparse(substitute(source)), "' exist",
			call. = FALSE)

# verify and adjust 'which.images' argument
	if (!is.null(which.images)) {
		if (length(which.images) == 2)
				which.images <- c(which.images, max(which.images))
		if (length(which.images) != 3)
			warning("'which.images' must be an integer vector of length 2 or 3",
			call. = FALSE)
		if (which.images[3] != max(which.images))
			warning("the third value in 'which.images' should be the largest",
			call. = FALSE)
	}

# assign 'method'
	if (method == "browser") {
		message("Some version of EBImage fails with the 'browser' option.")
		message("If this happens, change 'method' to 'raster'.")
		flush.console()
	if (method == "raster")
		opar <- par(ask = ask)
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
	pat1 <- "^[[:digit:]]{0,3}"
	pat2 <- "[abcdefghijklmnop][[:digit:]]+$"
	pat <- paste0(pat1, pat2)
	sel <- grepl(pat, field2, ignore.case = TRUE)

# assign value to imageType as "byWell" or "byFile" and complete message
	if (all(sel)) { # extract well and numeric optional prefix
		imageType <- "byWell"
		prefix <- sub(paste0("(^", pat1, ").*$"), "\\1", field2)
		well <- sub(paste0(pat1, "(.*$)"), "\\1", field2)
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
		ffsplit <- split(ff, field2)
	else if (imageType == "byFile")
		ffsplit <- split(ff, filename)
	else
		stop("CAN'T HAPPEN! Unexpected value for 'imageType'")

# read all images as a list and coerce to grayscale with a warning
	message("Reading images...", appendLF = FALSE)
	flush.console()
	img <- lapply(ffsplit, function(f) suppressWarnings(readImage(f)))
	message("passed")
	flush.console()
	if (any(sapply(img, colorMode) != 0)) {
		warning("images have been converted to grayscale by uniform RGB averaging",
			call. = FALSE)
		img <- lapply(img, channel, "gray")
	}

# perform check for expected numbers of images
	if (!is.null(which.images)) {
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
		message("Found ", nspl, " group", ifelse(nspl == 1, "", "s"), " of image pairs")
		flush.console()
	}

# report on and process each group of files in turn
	if (is.null(which.images)) { # simply show images in each file
		for (IDX in seq_along(ffsplit)) {
			nimg <- numberOfFrames(img[[IDX]])
			message(sprintf("%3d: %2d in %s", IDX, nimg, names(ffsplit)[IDX]))
			if (method != "none") {
				x <- normalize(img[[IDX]], separate = separate)
				opar <- par(ask = ask)
				display(x, all = TRUE, method = method, title = names(ffsplit)[IDX])
				par(opar)
			}
		}
	}
	else { # show DNA/target pairs
		for (IDX in seq_along(ffsplit)) {
			nimg <- numberOfFrames(dnaImages[[IDX]])
			message(sprintf("%3d: %2d in %s", IDX, nimg, names(ffsplit)[IDX]))
			if (method != "none") {
				myDna <- normalize(dnaImages[[IDX]], separate = separate)
				myTgt <- normalize(tgtImages[[IDX]], separate = separate)
				i <- c(rbind(seq_len(nimg), nimg + seq_len(nimg)))
#				img <- combine(myDna, myTgt)[,,i]
				img <- abind(myDna, myTgt, along = 3)[,,i]
				opar <- par(ask = ask)
				display(tile(img, nx = nx), all = TRUE, method = method,
					title = names(ffsplit)[IDX])
				par(opar)
			}
		}
	}
	invisible(ff) # return list of file names
}

