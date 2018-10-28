#' Extract Fluorescent Intensities from Paired Microscopic Images
#'
#' Extract the mean fluorescence intensity of DNA and a second fluorescent
#' target image for individual cells from paired images.
#'
#' @param source A character vector identifying a directory or directories
#'   with either multilayer tiff image files \emph{or} subdirectories identified
#'   by well with separate, paired images within each subdirectory for each
#'   well \emph{or} a character vector of image files organized as described.
#' @param type A character string identifying the type of image files to parse
#'   ("tif", "tiff", "jpeg", "jpg" or "png".)
#' @param which.images An integer of length 2 or 3. The first two numbers indicate
#'   the position of the DNA image and the target image. The optional
#'   third number specifies the total number of images for each field. A value of
#'   c(1, 2) indicates DNA first and target second. A value of c(2, 1) indicates
#'   target first and DNA image second. A value of c(1, 2, 4)
#'   indicates DNA image, target image, and two additional (ignored) images
#'   such as a second fluorescent color and phase contrast image for each field.
#' @param pattern Optional \code{grep} pattern as a character string used
#'   by \code{list.files()} to select image files.
#' @param args.nucMask A list of arguments passed to \code{nucMask()}.
#' @param args.trimMask A list of arguments passed to \code{trimMask()} such as
#'   \code{cutoff} or \code{k}. If this value is \code{NA}, no trimming is performed.
#' @param cellMask.flag If \code{TRUE}, the default
#'   nuclear mask will be used to generate a mask with \code{cellMask()}. This
#'   "cellular" mask will be used to measure fluorescence in the target image.
#' @param equalize \emph{If the fluorescent target images have more background
#'   pixels than foreground pixels} and if the background varies significantly
#'   from image to image, this can be set to \code{TRUE} in order to equalized
#'   the fluorescent images by subtracting the median value after applying a
#'   median filter and gaussian blur using the function \code{bnormalize()}.
#'
#' @details
#'
#' This is the core function to read and parse image data in a suite of tools
#' implemented with \code{\link{EBImage}} intended to determine
#' viral titers from sets of fluorescent micrographs. Typically, the images are
#' acquired in pairs where the first of each pair
#' is a DNA image and the second a fluorescent image of the viral target.
#' Because individual cells are identified by the nuclear stain, it
#' \emph{may} be beneficial to collect overexposed DNA images.
#'
#' These tools were developed to process fluorescent virus titers
#' performed in multi-well plates and is designed to parse images
#' collected at different multiplicities of infection or moi. The
#' moi can be expressed as virions (VP) per cell \emph{or} infectious
#' units (IU) per cell \emph{or} a unit of volume (ml, ul, nl) per cell.
#' These details are added to output of this function with the
#' \code{mergePdata()} function. The default order is for the nuclear
#' (typically DAPI) image to precede the viral antigen image, but this
#' sequence can be adjusted with the \code{which.images} argument.
#'
#' Images associated with each moi can be individual files in a
#' single directory where each directory is named for the well such as
#' \code{A1}, \code{A2}, etc. and the files within are identified as
#' \code{A1/file001.tif}, \code{A1/file002.tif}, etc. The well identifier
#' can be in upper or lower case and can contain leading zeros such as
#' \code{c0003/file12.tif}.
#'
#' Alternatively, each group of images associated with a given moi can be
#' a multi-layered tiff file where the sequence of images in the file is
#' specified by the argument \code{which.images}.
#'
#' A unique ID for each image can be created from a combination of
#' \code{frame} and either \code{well} or \code{file}. This can be useful
#' if it is necessary to determine a separate cutoff value for each target.
#' \preformatted{
#'   df$uid <- with(df, interaction(well, frame))
#'   df$uid <- with(df, interaction(file, frame))
#' }
#'
#' If the images have background values that vary from image to image or
#' have significant noise and if the fluorescent images have \emph{more}
#' background pixels than foreground pixels, then the argument
#' \code{equalize} can be set to \code{TRUE} to smooth the images by
#' sequentially modifying the values in each target image with a
#' median filter of radius 2, a Gaussian blur of radius 2, subtracting
#' the median value for each image and adding an offset of 0.05.
#' \emph{This cannot be used for images that have a large fraction
#' of positive cells, such as the example data set.}
#'
#' @return
#'
#' A data.frame of processed image data. \strong{All} data.frames will
#' have the following variables:
#' \describe{
#'   \item{\code{frame}}{Image sequence within each level of well
#'     or file as a factor (1, 2, 3, ...)}
#'   \item{\code{xm, ym}}{Center of mass (in pixels) for nucleus.}
#'   \item{\code{area}}{Area of mask used to calculate target mfi.}
#'   \item{\code{dna}}{Mean fluorescence intensity for DNA stain,
#'     typically not meaningful with over-exposed images.}
#'   \item{\code{mfi}}{Mean fluorescence intensity for the target,
#'     measured with selected mask.}
#' }
#' Results from data organized by \strong{well} will also include:
#' \describe{
#'   \item{\code{well}}{Harmonized well identifier from the \code{well.info()} function.}
#'   \item{\code{row}}{Row identifier ("A", "B", "C", etc.) as a factor.}
#'   \item{\code{column}}{Column number as a factor.}
#' }
#' Results from data organized as \strong{stacks} (multi-layered
#' tiff files) will include:
#' \describe{
#'   \item{\code{file}}{The file name as a factor.}
#' }
#'
#' @examples
#' # Note that execution of these examples can be rather slow...
#'   path.by.folder <- system.file("extdata", "by_folder", package = "virustiter")
#'   df.by.folder <- parseImages(path.by.folder)
#'   head(df.by.folder)
#'
#'   path.by.stack <- system.file("extdata", "by_stack", package = "virustiter")
#'   df.by.stack <- parseImages(path.by.stack)
#'   head(df.by.stack)
#'
#' # plots
#'   opar <- par(mfrow = c(1,2))
#'   plot(log(mfi) ~ area, df.by.folder, main = "By Folder", las = 1)
#'   plot(log(mfi) ~ area, df.by.stack, main = "By Stack", las = 1)
#'   par(opar)
#'
#' @import EBImage
#'
#' @export
#'
parseImages <- function(source, type = "tiff", which.images = c(1, 2, 2),
	pattern = NULL, args.nucMask = NULL, args.trimMask = NULL,
	cellMask.flag = FALSE, equalize = FALSE)
{
# requires EBImage, ensure appropriate values for parameters
	if (!require(EBImage))
		stop("The 'EBImage' package must be installed with biocLite")

# are all files or directories found in 'source' argument legitimate?
	if (!all(file.exists(source)))
		stop("not all files named in ", deparse(substitute(source)), " exist")

	if (all(file.info(source)$isdir))
		ff <- list.images(path = source, type = type, pattern = pattern)
	else if (all(!file.info(source)$isdir))
		ff <- source
	else
		stop("unable to use files/source in ", deparse(substitute(source)))

# check on arguments
	if (length(which.images) == 2)
			which.images <- c(which.images, max(which.images))
	if (length(which.images) != 3)
		stop("'which.images' must be an integer vector of length 2 or 3")
	if (which.images[3] != max(which.images))
		stop("the third value in 'which.images' must be the largest")

# Evaluate and read image data 
	message("Evaluating images...", appendLF = FALSE)

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
		stop("\nunable to use mixture of image files in ", source, '"')

# split image paths into related groups (by well or by file)
	if (imageType == "byWell")
		ffsplit <- split(ff, well)
	else if (imageType == "byStack")
		ffsplit <- split(ff, filename)
	else
		stop("\nunexpected value for 'imageType'")

# read all images as a list and coerce to grayscale with a warning
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
		stop("\nThe number of images in ", names(img)[bad],
			" is not a multiple of ", n_field)

# extract dna images
	dnaImages <- lapply(img, function(x, first = n_dna, by = n_field) {
		N <- dim(x)[3]
		if (N <= 1) stop("\nUNEXPECTED image dimension in dnaImages") 
		x[,,seq(first, N, by), drop = FALSE]})

# extract mfi images
	mfiImages <- lapply(img, function(x, first = n_mfi, by = n_field) {
		N <- dim(x)[3]
		if (N <= 1) stop("\nUNEXPECTED image dimension in mfiImages") 
		x[,,seq(first, N, by), drop = FALSE]})
	message("done")	# evaluating images

# option to smooth and equalize mfi images
	if (equalize == TRUE) {
		message("Equalizing target images...", appendLF = FALSE)
		mfiImages <- lapply(mfiImages, bnormalize)
		message("done")
	}

# initialize variable to collect results, initialize progress bar
	nff <- length(ffsplit)
	ret <- rep(list(NULL), nff)
	showProgress <- ifelse(nff > 2, TRUE, FALSE)

	message("Processing images.")
	if (showProgress)
		pb <- txtProgressBar(min = 1, max = nff, style = 3)

################################################################################
#
# process each group of files in turn
#
	for (IDX in seq_along(ffsplit)) {
		if (showProgress)
			setTxtProgressBar(pb, IDX)

	# create mask with any additional arguments in args.nucMask
		arg.list <- as.list(args(nucMask))
		arg.list <- arg.list[names(arg.list) != ""] # drop NULL values
		arg.list$dna <- dnaImages[[IDX]]
		nms <- names(args.nucMask)
		nms <- nms[nms %in% names(arg.list)] # find replacements
		sel <- names(arg.list) %in% nms
		arg.list <- c(arg.list[!sel], args.nucMask[nms])
		nmask <- do.call(nucMask, arg.list)

	# remove small and large nuclei with arguments in args.trimMask
		if (is.null(args.trimMask) || (!is.na(args.trimMask))) {
			arg.list <- as.list(args(trimMask))
			arg.list <- arg.list[names(arg.list) != ""]
			arg.list$mask <- nmask
			nms <- names(args.trimMask)
			nms <- nms[nms %in% names(arg.list)]
			sel <- names(arg.list) %in% nms
			arg.list <- c(arg.list[!sel], args.trimMask[nms])
			nmask <- do.call(trimMask, arg.list)
		}

	# expand mask to use estimated cell
		if (cellMask.flag == TRUE)
			cmask <- cellMask(nmask)
		else
			cmask <- nmask

	# ensure that images have three dimensions
		myDna <- dnaImages[[IDX]]
		myMfi <- mfiImages[[IDX]]
		dm <- dim(nmask)
		if (length(dm) == 2)
			 dim(myDna) <- dim(myMfi) <- dim(cmask) <- dim(nmask) <- c(dm, 1)
		nframes <- dim(nmask)[3]

	# measure and assemble data for the group indexed by 'IDX'
		res <- data.frame()
		for (i in seq_len(nframes)) {
			area <- computeFeatures.shape(cmask[,,i])[,1]
			XY <- computeFeatures.moment(nmask[,,i])[,1:2]
			dna <- computeFeatures.basic(nmask[,,i], myDna[,,i])[,1]
			mfi <- computeFeatures.basic(cmask[,,i], myMfi[,,i])[,1]
			if (imageType == "byWell") {
				ww <- names(ffsplit)[IDX]
				res <- rbind(res, data.frame(well = well.info(ww)$well,
					row = well.info(ww)$row,
					column = well.info(ww)$column, frame = i,
					xm = XY[,1], ym = XY[,2], area, dna, mfi))
			}
			else # imageType == "byStack"
				res <- rbind(res, data.frame(file = names(ffsplit)[IDX],
					frame = i, xm = XY[,1], ym = XY[,2],
					area, dna, mfi))
		}
		res$frame <- factor(res$frame, levels = sort(unique(res$frame)))
		rownames(res) <- NULL
		ret[[IDX]] <- res
	}
#
# done with processing each group of files
#
################################################################################

# compile and return collected data
	ans <- do.call(rbind, ret)
	rownames(ans) <- NULL
	if (showProgress)
		close(pb)
	message("Done")
	return(ans)
}
