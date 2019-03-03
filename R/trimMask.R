#' Remove Objects by Size, Proximity to Edge and Eccentricity
#' 
#' Remove large and small objects from an integer \code{Image} mask,
#' remove objects that are near the edge of the image, and remove
#' eccentric objects.
#' 
#' @param mask Object mask or list of masks with connected pixels having
#'   the same integer value.
#' @param cutoff Optional integer value of length 2 specifying the lower 
#'   and upper limits for the area in pixels. A value of \code{FALSE} or
#'   \code{NA} prevents size exclusion from occurring. A value of
#'   \code{NULL} makes use of the multiplier parameter \code{k} to determine 
#'   the cutoff limits. Either of the two values in \code{cutoff} can be
#'   specified as \code{NA} to use the multiplier parameter for that position.
#' @param k Numeric value of length 2 specifying the lower and upper
#'   multiplier to determine the cutoff from the \code{\link{median}} and
#'   \code{\link{mad}} of the area if \code{cutoff} is \code{NULL}.
#' @param border Objects within this many pixels of the edge will be excluded.
#' @param brush After coercion to the nearest odd integer, values > 0
#'   \code{\link[EBImage:morphology]{dilate}} the mask whereas values < 0
#'   \code{\link[EBImage:morphology]{erode}} the mask. 
#' @param ecc.max Exclude objects with elliptical eccentricity greater than
#'   this value.
#' 
#' @details
#' 
#' For each non-\code{NA} value in \code{cutoff}, objects smaller than
#' \code{cutoff[1]} and larger than \code{cutoff[2]} will be removed. 
#' Otherwise, objects smaller than \code{median(area) - k[1]*mad(area)} and
#' larger than \code{median(area)+ k[2]*mad(area)} will be removed. Objects
#' that are within \code{border} pixels of the edge will be removed. Objects
#' that have eccentricity greater than \code{ecc.max} will be removed. A circle
#' has eccentricity of 0 and a straight line has eccentricity of 1. Note 
#' that EBImage provides a poor approximation of this value. The final mask
#' will be dilated or eroded if \code{brush} is non-zero using a disc-shaped
#' brush with size determined by the adjusted value of \code{brush}.
#' 
#' @return
#'
#' Object mask or list of masks with objects removed and re-enumerated.
#'
#' @examples
#'
#'   x <- readImage(system.file("extdata", "by_folder/a1/file001.tif", package = "virustiter"))
#'   xm <- nucMask(x)
#'   xm2 <- trimMask(xm, cutoff = c(200, 400))
#'   xm3 <- trimMask(xm, border = 24)
#'   xm4 <- trimMask(xm, ecc.max = 0.75)
#'   plot(colorLabels(combine(xm, xm2, xm3, xm4)), all = TRUE)
#'   sapply(list(xm, xm2, xm3, xm4), apply, 3, max) # how many remain?
#'
#' @import EBImage
#' 
#' @export
#'
trimMask <- function(mask, cutoff = FALSE, k = c(1.5, 3), border = 0, brush = 0,
	ecc.max = 1)
{
	require(EBImage)
	if (missing(mask)) {
		usage <- c("trimMask argument hints:",
			'  k = c(3,5) to drop objects < 3x mad(area) and > 5x mad(area)',
			'  cutoff = NULL will  use values in k to trim',
      '  cutoff = c(100, Inf) to drop objects < 100 pixels',
      '  brush = -5 to erode mask with disc of radius 5 pixels',
      '  brush = 5 to dilate mask with disc of radius 5 pixels',
			'  border = 2 to drop objects within 2 pixels of image border',
      '  ecc.max = 0.75 to drop objects with eccentricity > 0.75')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}


# process function
	.proc <- function(mask, cutoff, k, border, brush, ecc.max)
	{
	# ensure three dimensions are present
		dm <- dim(mask)
		if (length(dm) == 2) dim(mask) <- c(dm, 1)
		nframes <- dim(mask)[3]
	# trim by area
		if (!identical(cutoff, FALSE) && !identical(cutoff, NA)) {
			area <- lapply(seq_len(nframes), function(i) computeFeatures.shape(mask[,,i])[,1])
			xmed <- median(unlist(area))
			xmad <- mad(unlist(area))
			if (is.null(cutoff)) cutoff <- c(NA, NA)
			if (is.na(cutoff[1])) cutoff[1] <- xmed - k[1] * xmad
			if (is.na(cutoff[2])) cutoff[2] <- xmed + k[2] * xmad
			lower <- max(cutoff[1], min(unlist(area)))
			upper <- min(cutoff[2], max(unlist(area)))
			small <- lapply(area, function(z) which(z < lower))
			large <- lapply(area, function(z) which(z > upper))
			mask <- rmObjects(mask, small, reenumerate = FALSE)
			mask <- rmObjects(mask, large, reenumerate = FALSE)
			mask <- reenumerate(mask) # needed because reenumeration may not OCCUR for no large
		}
	# trim border objects
		if (border > 0) {
			sel <- edgeObjects(mask, border = border)
			mask <- rmObjects(mask, sel, reenumerate = TRUE)
		}
	# trim by eccentricity
		if (ecc.max != 1) {
			ecc <- lapply(seq_len(nframes),
				function(i) computeFeatures.moment(mask[,,i])[,"m.eccentricity"])
			sel <- lapply(ecc, function(v) which(v > ecc.max))
			mask <- rmObjects(mask, sel, reenumerate = TRUE)
		}
	# apply erosion or dilation
		brush <- as.integer(brush)
		if (brush != 0) {
			brush <- 2*(brush + ifelse(brush < 0, -1, 0))%/%2 + 1 # ensure odd number
			if (brush > 0) # dilating preserves values of mask
				mask <- dilate(mask, makeBrush(brush, "disc"))
			else if (brush < 0) { # but eroding creates a binary object!
				mult <- erode(mask, makeBrush(-brush, "disc"))
				mask <- mult * mask # converts to eroded integer mask 
			}
		}
	# return mask to original dimensions 
		dim(mask) <- dm
		return(mask)
	}

# dispatch function according to argument 'mask'
	if (is(mask, "Image"))
		ans <- .proc(mask, cutoff = cutoff, k = k, border = border,
				brush = brush, ecc.max = ecc.max)
	else if (all(sapply(mask, is, "Image")))
		ans <- lapply(mask, .proc, cutoff = cutoff, k = k, border = border,
				brush = brush, ecc.max = ecc.max)
	else
		stop("'mask' must be an Image or list of images")
	return(ans)
}
