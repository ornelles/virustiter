#' Remove Objects by Size, Proximity to Edge and Eccentricity
#' 
#' Remove large and small objects from an integer \code{Image} mask,
#' remove objects that are near the edge of the image, and remove
#' highly eccentric objects.
#' 
#' @param mask Object mask or list of masks with connected pixels having
#'   the same integer value.
#' @param cutoff Optional integer value of length 2 specifying the lower 
#'   and upper limits for the area  in pixels. Either can be specified
#'   as \code{NA} to use the multiplier parameter, \code{k}.
#' @param k Numeric value of length 2 specifying the lower and upper
#'   multiplier to determine the cutoff from the \code{mean} and \code{mad}
#'  of the area.
#' @param border Exclude objects within this many pixels from the edge.  
#' @param ecc.max Exclude objects with elliptical eccentricity greater than
#'   this value.
#' 
#' @details
#' 
#' For each non-\code{NA} value in \code{cutoff}, objects smaller than
#' \code{cutoff[1]} and larger than \code{cutoff[2]} will be removed. 
#' Otherwise, objects smaller than \code{mean(area) - k[1]*mad(area)} and
#' larger than \code{mean(area)+ k[2]*mad(area)} will be removed. Objects
#' that are within \code{border} pixels of the edge will be removed. Objects
#' that have eccentricity greater than \code{ecc.max} will be removed. A circle
#' has eccentricity of 0 and a straight line has eccentricity of 1. Note 
#' that EBImage provides a poor approximation of this value.
#' 
#' @return
#'
#' Object mask or list of masks with objects removed.
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
trimMask <- function(mask, cutoff = NULL, k = c(1.5, 3), border = 0, ecc.max = 1)
{
	require(EBImage)
# process function
	.proc <- function(mask, cutoff, k, border, ecc.max)
	{
	# ensure three dimensions are present
		dm <- dim(mask)
		if (length(dm) == 2) dim(mask) <- c(dm, 1)
		nframes <- dim(mask)[3]
	# trim by area
		area <- apply(mask, 3, function(v) computeFeatures.shape(v)[,1])
		if (is(area, "matrix")) area <- split(area, c(col(area)))
		xbar <- mean(unlist(area))
		xmad <- mad(unlist(area))
		if (is.null(cutoff)) cutoff <- c(NA, NA)
		if (is.na(cutoff[1])) cutoff[1] <- xbar - k[1] * xmad
		if (is.na(cutoff[2])) cutoff[2] <- xbar + k[2] * xmad
		lower <- max(cutoff[1], min(unlist(area)))
		upper <- min(cutoff[2], max(unlist(area)))
		small <- lapply(area, function(z) which(z < lower))
		large <- lapply(area, function(z) which(z > upper))
		mask <- rmObjects(mask, small, reenumerate = FALSE)
		mask <- rmObjects(mask, large)
	# trim border objects
		if (border > 0) {
			sel <- edgeObjects(mask, border = border)
			mask <- rmObjects(mask, sel)
		}
	# trim by eccentricity
		if (ecc.max != 1) {
			ecc <- apply(mask, 3, function(v) computeFeatures.moment(v)[,"m.eccentricity"])
			if (is(ecc, "matrix")) ecc <- split(ecc, c(col(ecc)))
			sel <- lapply(ecc, function(v) which(v > ecc.max))
			mask <- rmObjects(mask, sel)
		}
		dim(mask) <- dm
		return(mask)
	}

# dispatch function according to argument 'mask'
	if (is(mask, "Image"))
		ans <- .proc(mask, cutoff = cutoff, k = k, border = border, ecc.max = ecc.max)
	else if (all(sapply(mask, is, "Image")))
		ans <- lapply(mask, .proc, cutoff = cutoff, k = k, border = border, ecc.max = ecc.max)
	else
		stop("'mask' must be an Image or list of images")
	return(ans)
}
