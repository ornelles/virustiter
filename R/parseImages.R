#' Parse Paired Microscopic Images
#'
#' Identify nuclei by a DNA stain and measure the fluorescence intensity of
#' the DNA and a second second fluorescent target image from paired images.
#'
#' @param nuc A list of nuclear images. If the second argument \code{tgt} is
#'   \code{NULL}, the first argument \code{nuc} is treated as a list of length
#'   2 containing  nuclear and target images.
#' @param tgt A list of fluorescent images corresponding to the nuclear images
#'   in \code{nuc}. If this argument is \code{NULL}, \code{nuc} must be a
#'   list of both image types. 
#' @param args.nucMask A list of arguments passed to \code{\link{nucMask()}}.
#' @param args.trimMask A list of arguments passed to \code{\link{trimMask()}}
#'   such as\code{cutoff} or \code{k}. If this value is \code{NA}, no
#'   trimming is performed.
#' @param cellMask.flag If \code{TRUE}, the default
#'   nuclear mask will be used to generate a mask with \code{cellMask()}. This
#'   larger mask will be used to measure fluorescence in the target image.
#' @param equalize If the fluorescent target images have \emph{more background
#'   pixels than foreground pixels} and if the background varies significantly
#'   from image to image, this can be set to \code{TRUE} in order to equalized
#'   the fluorescent images by subtracting the median value after applying a
#'   median filter and gaussian blur using the function \code{bnormalize()}.
#'
#' @details
#'
#' This function identifies cells by a DNA stain and measures the
#' fluorescent intensity in both the DNA stain and paired fluorescent image.
#' This is part of a suite of tools implemented with \code{\link{EBImage}}
#' designed to determine viral titers from sets of fluorescent micrographs.
#' The first argument to this function can be the result of the function
#' \code{\link{getImages}}. Images provided to this function are pairs where
#' the first of each pair is a DNA image and the second is a fluorescent image
#' of the viral antigen. Because individual cells are identified by the nuclear
#' stain, it is often beneficial to collect overexposed DNA images.
#'
#' These tools were developed to process fluorescent virus titers
#' performed in multi-well plates and is designed to parse images
#' collected at different multiplicities of infection or moi. The
#' moi can be expressed as virions (VP) per cell \emph{or} infectious
#' units (IU) per cell \emph{or} a unit of volume (ml, ul, nl) per cell.
#' These details are added to output of this function with the
#' \code{mergePdata()} function. The required order is for the nuclear
#' (typically DAPI) image to precede the viral antigen image. This
#' sequence can be adjusted with \code{which.images} argument in
#' the function \code{\link{getImages}}.
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
#' if it is necessary to determine a separate cutoff value for each pair of
#' images.
#' \preformatted{
#'   df$uid <- with(df, interaction(well, frame))
#'   df$uid <- with(df, interaction(file, frame))
#' }
#'
#' If the images have background values that vary from image to image or
#' have significant noise and if the fluorescent images have \emph{more}
#' background pixels than foreground pixels, then the argument
#' \code{equalize} can be set to \code{TRUE}. This smooths the images by
#' sequentially modifying values in each target image with a
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
#'   \item{\code{area}}{Area of the mask (in pixels) used to calculate target mfi.}
#'   \item{\code{dna}}{Mean fluorescence intensity for DNA stain,
#'     typically not meaningful if the DNA image was over-exposed.}
#'   \item{\code{mfi}}{Mean fluorescence intensity for the target,
#'     measured with selected mask.}
#' }
#' Results from data organized by \strong{well} will also include:
#' \describe{
#'   \item{\code{well}}{Harmonized well identifier from the \code{well.info()} function.}
#'   \item{\code{row}}{Row identifier ("A", "B", "C", etc.) as a factor.}
#'   \item{\code{column}}{Column number as a factor.}
#' }
#' Results from data organized as multi-layered tiff files (or stacks)
#' will include:
#' \describe{
#'   \item{\code{file}}{The file name as a factor.}
#' }
#'
#' @examples
#' # Note that execution of these examples can be rather slow...
#'   path.by.folder <- system.file("extdata", "by_folder", package = "virustiter")
#'   images.by.folder <- getImages(path.by.folder)
#'   df.by.folder <- parseImages(images.by.folder)
#'   head(df.by.folder)
#'
#'   path.by.stack <- system.file("extdata", "by_stack", package = "virustiter")
#'   images.by.stack <- getImages(path.by.stack)
#'   df.by.stack <- parseImages(images.by.stack)
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
parseImages <- function(nuc, tgt = NULL, args.nucMask = NULL,
	args.trimMask = NULL, cellMask.flag = FALSE, equalize = FALSE)
{
# requires EBImage, ensure appropriate values for parameters
	if (!require(EBImage))
		stop("The 'EBImage' package must be installed with biocLite")

# check on arguments
	if (is.null(tgt)) {
		if (!is.list(nuc) || length(nuc) != 2)
			stop("'nuc' must be a list oflength two")
		dnaImages <- nuc[[1]]
		mfiImages <- nuc[[2]]
	}
	else if (is(nuc, "Image") & is(tgt, "Image")) {
		dnaImages <- list(nuc)
		mfiImages <- list(tgt)
	}
	else {
		if (!is.list(nuc) || !is.list(tgt))
			stop("'nuc' and 'tgt' must be Image objects or lists of Image objects")
		dnaImages <- nuc
		mfiImages <- tgt
	}

# determine imageType as "well" or "file"
	sel <- grepl("^[abcdefghijklmnop][[:digit:]]+$", names(dnaImages),
		ignore.case = TRUE)
	if (length(sel) > 0 && all(sel))
		imageType <- "byWell"
	else if (length(sel) > 0 && all(!sel))
		imageType <- "byFile"
	else {
		imageType <- "byFile"
		names(dnaImages) <- sprintf("image%04d", seq_along(dnaImages))
		Message("Unable to determine organization of imagesget")
	}
 
# option to smooth and equalize mfi images
	if (equalize == TRUE) {
		message("Equalizing target images...", appendLF = FALSE)
		mfiImages <- lapply(mfiImages, bnormalize)
		message("done")
	}

# initialize variable to collect results, initialize progress bar
	nff <- length(dnaImages)
	ret <- rep(list(NULL), nff)
	showProgress <- ifelse(nff > 2, TRUE, FALSE)

	message("Processing images by ", ifelse(imageType == "byWell", "well", "file"))
	if (showProgress)
		pb <- txtProgressBar(min = 1, max = nff, style = 3)

################################################################################
#
# process each group of files in turn
#
	for (IDX in seq_len(nff)) {
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
				ww <- names(dnaImages)[IDX]
				res <- rbind(res, data.frame(well = well.info(ww)$well,
					row = well.info(ww)$row,
					column = well.info(ww)$column, frame = i,
					xm = XY[,1], ym = XY[,2], area, dna, mfi))
			}
			else # imageType == "byStack"
				res <- rbind(res, data.frame(file = names(dnaImages)[IDX],
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
