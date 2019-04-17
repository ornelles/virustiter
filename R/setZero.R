#' Normalize Image to Given Zero Value
#'
#' Scale an image and crop at the given zero value.
#'
#' @param x A numeric vector or array \emph{or} list of such objects, where
#'   each object is typically a grayscale \code{Image} object.
#' @param zero Numeric value or list of values specifying the
#'   zero value pixel for each frame.
#' @param min.value Numeric value or list of values to be the
#'   new minimum value in the transformed image, default of 0.
#'
#' @details
#' 
#' Each frame of the argument will be linearly scaled by subtracting
#' \code{zero}. Values less than \code{min.value} will be set to
#' \code{min.value}. Values of \code{0, NULL, or NA} for \code{min.value}
#' will be treated as the value 0. 
#'
#' This function will probably replace \code{\link{bnormalize}}
#' soon. This function is typically applied after using
#' \code{\link{getZero}} to determine the likely zero value.
#'
#' @seealso \code{\link{bnormalize}}, \code{\link{getZero}}
#'
#' @return
#'
#' An object of the same structure as the argument \code{x} that has been
#' linearly scaled and cropped.
#'
#' @export
#'
setZero <- function(x, zero, min.value = 0)
{
# argument check
	if (missing(min.value) || is.na(min.value) || is.null(min.value))
		min.value <- 0
	if (min.value > 1 | min.value < 0)
		warning("'min.value' outside of the range [0,1]")

# working function
	.setZero <- function(x, zero, min.value) {
		x <- x - zero
		x[x < min.value] <- min.value
		return(x)
	}

# dispatch
	if (is(x, "list") & all(sapply(x, is.numeric)))
		ans <- Map(.setZero, x, zero, min.value)
	else if (is.numeric(x))
		ans <- .setZero(x, zero, min.value)
	else
		stop("require a numeric vector, numeric array or list of such objects")
	return(ans)
}
