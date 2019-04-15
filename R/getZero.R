#' Find Background (Zero Value) Pixel
#'
#' Find the most common non-zero pixel in the lower fraction
#' of each frame of an image.
#'
#' @param img An array object or list of array objects, where each
#'   is typically a grayscale \code{Image} object.
#' @param frac Number between 0 and 1 specifying the fraction of
#'   values to consider.
#'
#' @details
#' 
#' This function is meant to be applied to images that have more
#' background pixels than foreground pixels. The most common
#' value among the lower half (as specified by \code{frac})
#' of intensity values is sought as the true 'zero' value.
#' 
#' Intensity values greater than the minimum and less than
#' \code{frac * maximum} are processed by the \code{\link{density}}
#' function to determine the most common pixel. The minimum pixel
#' value is excluded to avoid skewing the density estimate if an
#' undersaturated image is used.
#'
#' @seealso \code{\link{bnormalize}}, \code{\link{setZero}},
#'   \code{\link{getBgnd}}, \code{\link{findBgnd}}
#'
#' @return
#'
#' A numeric vector or list of numeric vectors with the most
#' common non-zero pixel for each frame.
#'
#' @export
#'
getZero <- function(img, frac = 0.5)
{
# argument check
	stopifnot(frac <= 1, frac >= 0)

# working function
	.getZero <- function(img, frac) {
		dm <- dim(img)
		if (length(dm) == 2) dim(img) <- c(dm, 1)
		xmax <- apply(img, 3, function(v) frac*diff(range(v)))
		xmin <- apply(img, 3, min)
		sapply(seq_len(dim(img)[3]), function(i) {
			valid <- img[,,i] > xmin[i] & img[,,i] < xmax[i]
			d <- density(img[,,i][valid])
			d$x[which.max(d$y)]})
	}

# dispatch
	# dispatch
	if (is(img, "list") & all(sapply(img, is.array)))
		ans <- lapply(img, .getZero, frac = frac)
	else if (is(img, "array"))
		ans <- .getZero(img, frac = frac)
	else
		stop("require an array or list of array objects")
	return(ans)
}
