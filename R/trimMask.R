#' Remove Objects by Size
#' 
#' Remove objects from an integer Image mask with areas outside of lower
#' and upper cutoff values or by a lower and upper multiplier applied to
#' the \code{mad} of the area. 
#' 
#' @param mask Object mask with connected pixels having the same
#' integer value.
#' @param cutoff Optional integer value of length 2 specifying the lower 
#' and upper absolute cutoff values in pixels.
#' @param k Numeric value of length 2 specifying the lower and upper multiplier 
#' to determine the cutoff from the \code{mean} and \code{mad} of the area.
#' @param reenumerate Re-enumerate the objects in the trimmed mask.
#' 
#' @details
#' 
#' If \code{cutoff} is specified, objects smaller than \code{cutoff[1]} and 
#' larger than \code{cutoff[2]} will be removed. Otherwise, objects smaller 
#' than \code{mean(area) - k[1]*mad(area)} and larger than \code{mean(area)+ 
#' k[2]*mad(area)} will be removed. The mask will be reenumerated if that 
#' parameter is \code{TRUE}. 
#' 
#' @return
#'
#' Object mask with small and large objects removed.
#'
#' @examples
#'
#' x <- readImage(system.file("extdata", "by_folder/b2/file002.tif", package = "virustiter"))
#' xb <- normalize(gblur(x, 2))
#' xt <- thresh(xb)
#' xm <- bwlabel(xt)
#' xm2 <- trimMask(xm)
#' xm3 <- trimMask(xm, cutoff = c(1, 20))
#' xm4 <- trimMask(xm, cutoff = c(100, 250))
#' plot(combine(xm, xm2, xm3, xm4), all = TRUE)
#'
#' @import EBImage
#' 
#' @export

trimMask <- function(mask, cutoff = NULL, k = c(1.5, 3), reenumerate = TRUE)
{
	require(EBImage)
	dm <- dim(mask)
	if(length(dm) < 2 || length(dm) > 3)
		stop("'mask' must be a 2- or 3-dimension integer array")
	if (length(dm) == 2)
		dim(mask) <- c(dm, 1)
	nframes <- dim(mask)[3]
	area <- lapply(1:nframes, function(i) computeFeatures.shape(mask[,,i])[,1])
	xbar <- mean(unlist(area))
	xmad <- mad(unlist(area))
	if (is.null(cutoff)) {
		lower <- xbar - k[1] * xmad
		upper <- xbar + k[2] * xmad
	}
	else {
		lower <- max(cutoff[1], min(unlist(area)))
		upper <- min(cutoff[2], max(unlist(area)))
	}
	small <- lapply(area, function(z) which(z < lower))
	large <- lapply(area, function(z) which(z > upper))
	mask <- rmObjects(mask, small, reenumerate = FALSE)
	mask <- rmObjects(mask, large)
	dim(mask) <- dm
	if (reenumerate)
		mask <- reenumerate(mask)
	return(mask)
}
