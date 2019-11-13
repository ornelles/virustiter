#' Find Baseline (Zero) Pixel Value
#'
#' Find the most common non-zero pixel in each frame of an image.
#'
#' @param x A numeric vector or array \emph{or} list of such objects, where
#'   each object is typically a grayscale \code{Image} object.
#' @param frac Number between 0 and 1 specifying the fraction of
#'   values to include in the search for the most common value.
#'
#' @details
#' 
#' Identify the  most common pixel value in each frame of an image
#' to serve as the true 'zero' value.
#' 
#' Intensity values greater than the minimum and less than
#' \code{frac * maximum} are processed by the \code{\link{density}}
#' function to determine the most common pixel. The minimum pixel
#' value is excluded to avoid skewing the density estimate if an
#' undersaturated image is used. If the image has many bright pixels,
#' limit the search to the lower fraction of values by setting
#' \code{frac} to a value such as \code{0.5}. 
#'
#' @seealso \code{\link{bnormalize}}, \code{\link{setZero}},
#'   \code{\link{getBgnd}}, \code{\link{findBgnd}}
#'
#' @return
#'
#' A numeric vector or list of numeric vectors with the most
#' common non-zero pixel for each frame.
#'
#' @examples
#' set.seed(123)
#' z <- c(rnorm(200, 1, 0.1), rnorm(200, 2, 1))
#' plot(density(z))
#' abline(v = getZero(z), col = 2)
#'
#' @export
#'
getZero <- function(x, frac = 1)
{
# argument check
	stopifnot(frac <= 1, frac >= 0)

# working function
	.getZero <- function(x, frac) {
		dm <- dim(x)
		if (is.null(dm)) dim(x) <- c(length(x), 1, 1)
		else if (length(dm) == 2) dim(x) <- c(dm, 1)
		xmin <- apply(x, 3, min)
		xdif <- apply(x, 3, function(v) frac*diff(range(v)))
		xmax <- xmin + xdif
		sapply(seq_len(dim(x)[3]), function(i) {
			valid <- x[,,i] > xmin[i] & x[,,i] < xmax[i]
			if (sum(valid) == 0) { # diff(range) is 0
				warning("frame ", i, " has diff(range) of 0", call. = FALSE)
				return(median(x[,,i]))
			}
			else {
				d <- density(x[,,i][valid])
				return(d$x[which.max(d$y)])
			}
		})
	}

# dispatch
	if (is(x, "list") & all(sapply(x, is.numeric)))
		ans <- lapply(x, .getZero, frac = frac)
	else if (is.numeric(x))
		ans <- .getZero(x, frac = frac)
	else
		stop("require a numeric vector, numeric array or list of such objects")
	return(ans)
}
