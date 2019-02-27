#' Linear Spatial Translation
#' 
#' Replacement for \code{EBImage} function that accepts lists
#' 
#' @param x An \code{Image} object or array or list of the same.
#' @param adj A vector of 2 numbers denoting the translation vector
#'   in pixels or a list of the same. 
#' @param filter A character string indicating the interpolating
#'   sampling filter. Valid values are \code{"none"} or the default, 
#'   \code{"bilinear"}.
#' @param ... Arguments to be passed to \code{affine} in the
#'   \code{EBImage} package, such as \code{filter}, \code{output.dim},
#'   \code{bg.col} or \code{antialias}.
#' 
#' @details
#' 
#' This function is a replacement for the EBImage code of the same name. 
#' 
#' @import EBImage
#' 
#' @export
#'
translate <- function (x, adj, filter = c("bilinear", "none"), ...) 
{
	filter <- match.arg(filter)
# ensure adj is a list
	if(!is.list(adj))
		adj <- list(adj)
# ensure image has 3 dimensions
	dm <- dim(x)
	if (length(dm) == 2)
		dim(x) <- c(dm, 1)
# check on compatibility of image and translations
	if(length(adj) != dim(x)[3])
		stop("'adj' is not the correct length for 'x'")
# process as three dimensional array
	m <- rbind(c(1, 0), c(0, 1), unlist(adj))
	dim(m) <- c(3L, 2L, dim(x)[3])
	ans <- Image(NA, dim = dim(x))
	for (i in seq_len(dim(x)[3]))
		ans[,,i] <- affine(x = x[,,i], m = m[,,i], filter = filter, ...)
	dim(ans) <- dm
	return(ans)
}
