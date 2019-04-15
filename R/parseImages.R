#' Parse Paired Microscopic Images
#'
#' Identify nuclei by a DNA stain and measure the fluorescence intensity of
#' the DNA and a second second fluorescent target image from paired images.
#'
#' @param nuc An \code{Image} object or list of such objects representing
#'   nuclear images. If the second argument (\code{tgt}) is \code{NULL},
#'   this must be a list of length 2 containing nuclear and target images.
#' @param tgt An \code{Image} object or list of fluorescent images corresponding
#'   to the nuclear images in \code{nuc}. If this argument is \code{NULL},
#'   \code{nuc} is treated as a list of both image types.
#' @param nMask An optional \code{Image} object or a list of
#'   these objects containing an integer \code{Image} mask identifying nuclei.
#'   If this value is \code{NULL}, the nuclear mask will be determined by
#'   \code{\link{nucMask}} with the arguments provided in \code{args.nMask}.
#' @param args.nMask A list of arguments passed to \code{\link{nucMask}}.
#'   This argument is ignored if a \code{nMask} is provided.
#' @param args.trimMask A list of arguments passed to \code{\link{trimMask}}.
#'   This argument is ignored if a \code{nMask} is provided. Otherwise, if
#'   \code{NULL}, \code{\link{trimMask}} will be called with default parameters.
#'   If \code{FALSE}, no trimming will be performed.
#' @param cMask An optional \code{Image} object or a list of 
#'   these objects containing an integer \code{Image} mask identifying
#'   regions of the target image in which to define cell boundaries. If 
#'   \code{TRUE}, the nuclear mask will be used to generate an expanded mask
#'   with \code{\link{cellMask}}. This larger mask will be used to measure
#'   fluorescence intensity in the target image. If \code{FALSE} (default),
#'   the nuclear mask will be used to measure fluorescence intensities. 
#' @param equalize If the background varies significantly among the fluorescent
#'   target images (or among the frames of target images), when \code{TRUE},
#'   this option adjusts the fluorescent target images to have a common
#'   background value with the default options in \code{\link{bnormalize}}
#'   and using the entire range of target images for the expected range.
#' @param simplify Return a single \code{data.frame} of results if \code{TRUE},
#'   otherwise return a list of \code{data.frames} for each member of the list.
#'
#' @details
#'
#' This function identifies cells by a DNA stain and measures the
#' fluorescent intensity in both the DNA image and the paired fluorescent image.
#' This is part of a suite of tools implemented with \code{\link{EBImage}}
#' designed to determine viral titers from sets of fluorescent micrographs.
#' The first argument to this function can be the result of the function
#' \code{\link{getImages}}. Images provided to this function are paired where
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
#' \code{\link{mergePdata}} function. Images must be ordered with the nuclear
#' (typically DAPI) image before the viral antigen image. Note that this
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
#' if it is necessary to determine a separate background value for each
#' pair of images.
#' \preformatted{
#'   df$uid <- with(df, interaction(well, frame))
#'   df$uid <- with(df, interaction(file, frame))
#' }
#'
#' @return
#'
#' A data.frame (or list of data.frames) containing  processed image
#' data. \strong{Each} data.frame will have the following variables:
#' \describe{
#'   \item{\code{frame}}{Image sequence within each level of well
#'     or file as an \code{integer} (1, 2, 3, ...)}
#'   \item{\code{xm, ym}}{Center of mass (in pixels) for nucleus.}
#'   \item{\code{area}}{Area of the mask (in pixels) used to calculate target mfi.}
#'   \item{\code{dna}}{Mean fluorescence intensity for DNA stain,
#'     typically not meaningful if the DNA image was over-exposed.}
#'   \item{\code{mfi}}{Mean fluorescence intensity for the target,
#'     measured with selected mask.}
#' }
#' Results from data organized by \strong{well} will also include:
#' \describe{
#'   \item{\code{well}}{Harmonized well identifier from the \code{\link{well.info}} function.}
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
parseImages <- function(nuc, tgt = NULL, nMask = NULL, cMask = FALSE,
	args.nMask = NULL, args.trimMask = NULL, args.cMask = NULL,
	equalize = FALSE, simplify = TRUE)
{
# requires EBImage, ensure appropriate values for parameters
	if (!require(EBImage))
		stop("The 'EBImage' package must be installed with biocLite")

# check first two arguments, extract images, ensure that they are lists
	if (is.null(tgt)) { # first argument is a list of length 2
		if (!is.list(nuc) || length(nuc) != 2)
			stop("'nuc' must be a list of length two")
		nucImages <- nuc[[1]]
		tgtImages <- nuc[[2]]
	}
	else if (is(nuc, "Image") & is(tgt, "Image")) { # both are images
		nucImages <- list(nuc)
		tgtImages <- list(tgt)
	}
	else { # both must be lists
		if (!is.list(nuc) || !is.list(tgt))
			stop("'nuc' and 'tgt' must be Image objects or lists of Image objects")
		nucImages <- nuc
		tgtImages <- tgt
	}

# determine imageType as "well" or "file"
	sel <- grepl("^[abcdefghijklmnop][[:digit:]]+$", names(nucImages),
		ignore.case = TRUE)
	if (length(sel) > 0 && all(sel))
		imageType <- "byWell"
	else if (length(sel) > 0 && all(!sel))
		imageType <- "byFile"
	else {
		imageType <- "byFile"
		names(nucImages) <- sprintf("image%04d", seq_along(nucImages))
		message("Unable to determine organization of images, using 'file'")
	}
 
# optionally equalize tgt images using the range of all tgt values
	if (equalize == TRUE) {
		message("Equalizing target images...", appendLF = FALSE)
		tgtRange <- range(sapply(tgt, range))
		tgtImages <- lapply(tgtImages, bnormalize, inputRange = tgtRange)
		message("done"); flush.console()
	}

# process nMask with additional arguments in args.nMask
	if (is.null(nMask)) { 
		message("Creating nuclear masks...", appendLF = FALSE)
		arg.list <- formals("nucMask")
		arg.list$dna <- nucImages
		nms <- names(args.nMask)
		nms <- nms[nms %in% names(arg.list)] # find replacements
		sel <- names(arg.list) %in% nms
		arg.list <- c(arg.list[!sel], args.nMask[nms])
		nmask <- do.call("nucMask", arg.list)
		message("done"); flush.console()
	}
	else {
		message("Using nuclear masks in '", deparse(substitute(nMask)), "'")
		nmask <- nMask
	}

# remove small and large nuclei with arguments in args.trimMask
	if (is.null(nMask) && (is.null(args.trimMask) ||
			(is.logical(args.trimMask) && args.trimMask == TRUE))) {
		message("Trimming nuclear masks...", appendLF = FALSE)
		arg.list <- formals("trimMask")
		arg.list$mask <- nmask
		nms <- names(args.trimMask)
		nms <- nms[nms %in% names(arg.list)]
		sel <- names(arg.list) %in% nms
		arg.list <- c(arg.list[!sel], args.trimMask[nms])
		nmask <- do.call("trimMask", arg.list)
		message("done"); flush.console()
	}

# process cMask
	if (is.logical(cMask) && cMask == TRUE) {
		message("Creating cell masks...", appendLF = FALSE)
		arg.list <- formals("cellMask")
		nms <- names(args.cMask)
		nms <- nms[nms %in% names(arg.list)] # find replacements
		sel <- names(arg.list) %in% nms
		arg.list <- c(arg.list[!sel], args.cMask[nms])
		cmask <- Map(cellMask, nmask, MoreArgs = arg.list)
		message("done"); flush.console()
	}
	else if(is.logical(cMask) && cMask == FALSE)
		cmask <- nmask
	else if (is(cMask, "list") || is(cMask, "Image")) {
		message("Using cell masks in '", deparse(substitute(cMask)), "'")
		cmask <- cMask
	}
	else
		stop("unable to use value in '",deparse(substitute(cMask)), "'")

# ensure that nmask and cmask are lists
	if (!is.list(cmask)) cmask <- list(cmask)
	if (!is.list(nmask)) nmask <- list(nmask)

# check that images and masks are compatible
	dm.nuc <- unname(sapply(nucImages, dim))
	dm.tgt <- unname(sapply(tgtImages, dim))
	dm.nm <- unname(sapply(nmask, dim))
	dm.cm <- unname(sapply(cmask, dim))

	if (!identical(dm.nuc, dm.tgt)) stop("nuc and tgt are mismatched")
	if (!identical(dm.nuc, dm.nm)) stop("nuc and nMask are mismatched")
	if (!identical(dm.nuc, dm.cm)) stop("nuc and cMask are mismatched")
	if (!identical(dm.nm, dm.cm)) stop("uh oh...this can't happen")

# initialize variable to collect results
	nImages <- length(nucImages)
	ans <- rep(list(NULL), nImages)

# initialize progress bar
	if (imageType == "byWell")
		message("Processing images by 'well'")
	else
		message("Processing images by 'file'")
	flush.console()
	showProgress <- ifelse(nImages > 1, TRUE, FALSE)
	if (showProgress)
		pb <- txtProgressBar(min = 1, max = nImages, style = 3)

################################################################################
#
# process each group of images - FAILING if one in set is two dimensional...
#
	for (idx in seq_len(nImages)) {
		if (showProgress) setTxtProgressBar(pb, idx)

	# extract appropriate element from lists
		myNuc <- nucImages[[idx]]; myTgt <- tgtImages[[idx]]
		myNmask <- nmask[[idx]]; myCmask <- cmask[[idx]]

	# ensure that images have three dimensions
		dm <- dim(myNuc)
		if (length(dm) == 2)
			 dim(myNuc) <- dim(myTgt) <- dim(myCmask) <- dim(myNmask) <- c(dm, 1)
		nframes <- dim(myNmask)[3]

	# measure and assemble data for the group indexed by 'idx'
		res <- data.frame()
		for (i in seq_len(nframes)) {
			area <- computeFeatures.shape(myCmask[,,i])[,1]
			XY <- computeFeatures.moment(myNmask[,,i])[,1:2]
			xm <- if(is.null(XY)) NULL else XY[,1]
			ym <- if(is.null(XY)) NULL else XY[,2]
			dna <- computeFeatures.basic(myNmask[,,i], myNuc[,,i])[,1]
			mfi <- computeFeatures.basic(myCmask[,,i], myTgt[,,i])[,1]
			if (imageType == "byWell") {
				ww <- names(nucImages)[idx]
				res <- rbind(res, data.frame(well = well.info(ww)$well,
					row = well.info(ww)$row, column = well.info(ww)$column,
					frame = i, xm = xm, ym = ym, area, dna, mfi))
			}
			else # imageType == "byFile"
				res <- rbind(res, data.frame(file = names(nucImages)[idx],
					frame = i, xm = xm, ym = ym,
					area, dna, mfi))
		}
		rownames(res) <- NULL
		ans[[idx]] <- res
	}
	names(ans) <- names(nucImages)
#
# done with processing each group of images
#
################################################################################

# compile and return collected data
	if (simplify == TRUE) {
		ans <- do.call("rbind", ans)
		rownames(ans) <- NULL
	}
	if (showProgress) close(pb)
	message("Done")
	return(ans)
}
