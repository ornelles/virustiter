#' Normalize Image to Given Zero Value
#'
#' Linearly scale an image to re-zero at the given zero value.
#'
#' @param img An array object or list of array objects, where each
#'   is typically a grayscale \code{Image} object.
#' @param zero Number between 0 and 1 for the new zero value.
#' @param bit.depth Integer between 1 and 64 specifying the
#'   bit-depth of the image, typically 8 or 16.
#' @param nonzero Logical value, if \code{TRUE}, the minimum value will
#'  adjusted to \code{2^-bit.depth}, otherwise the minimum value
#'  will be set to 0.
#'
#' @details
#' 
#' Each frame of the argument will be linearly scaled to have a
#' minimum at \code{zero} by subtraction. Negative values will
#' be clipped to 0 if \code{nonzero} is \code{FALSE} otherwise the
#' minimum value will be set to \code{2^-bit.depth}. 
#'
#' This function will probably replace \code{\link{bnormalize}}
#' soon. This function is typically applied after using
#' \code{\link{getZero}} to determine the likely zero value.
#'
#' @seealso \code{\link{bnormalize}}, \code{\link{getZero}}
#'
#' @return
#'
#' An array or list of the same size that has been linearly scaled.
#'
#' @export
#'
setZero <- function(img, zero, bit.depth = 16, nonzero = TRUE)
{
# argument check
	stopifnot(bit.depth >= 1, bit.depth <= 64)

# working function
	.setZero <- function(img, zero, bit.depth, nonzero) {
		img <- img - rep(zero, each = prod(dim(img)[1:2]))
		min.value <- if(nonzero) 1/2^as.integer(bit.depth) else 0
		img[img < zero] <- min.value
		return(img)
	}

# dispatch
	if (is(img, "list") & all(sapply(img, is.array)))
		ans <- Map(.setZero, img, zero, bit.depth = bit.depth, nonzero = nonzero)
	else if (is(img, "array"))
		ans <- .setZero(img, zero, bit.depth, nonzero)
	else
		stop("require an array or list of array objects")
	return(ans)
}
