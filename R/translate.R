#' Linear Spatial Translation
#' 
#' Replacement for \code{EBImage} function that accepts lists
#' 
#' @param x An \code{Image} object or array or list of the same.
#' @param v A vector of 2 numbers denoting the translation vector
#'   in pixels or a list of the same. 
#' @param filter A character string indicating the interpolating
#'   sampling filter. Valid values are \code{"none"} or the default, 
#'   which has been changed to \code{"bilinear"}.
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
#' @examples
#' 
#' x <- readImage(system.file("images", "sample-color.png", package = "EBImage"))
#' y <- untile(x, c(4, 4), lwd = 0)
#' plot(y, all = TRUE, nx = 4)
#' v <- sample(-20:20, 4*4*2, replace = TRUE)
#' v <- split(v, gl(4*4, 2)) # 16 pairs of translations
#' plot(translate(y, v), all = TRUE, nx = 4)
#' 
#' @export
#'
translate <- function (x, v, filter = c("bilinear", "none"), ...) 
{
	filter <- match.arg(filter)
# determine number of rendering frames
	nf <- numberOfFrames(x, type = "render")
# check for match between x and v
	if (nf == 1) {
		if (is.numeric(v))
			v <- list(v)
		v <- v[1]
	}
# adjust v if nf > 1
	if (nf > 1) {
		if (!is.list(v))
			stop("'v' must be a list of length ", nf, " for 'x'")
		else if (length(v) == 1) 
			v <- rep(v, nf)	
		else if (length(v) != nf){
			v <- rep(v, nf)[1:nf]
			warning("'v' replicated an uneven number of times")
		}
	}
# process each frame
	m <- rbind(c(1, 0), c(0, 1), unlist(v)) 
	dim(m) <- c(3L, 2L, nf) # 3 x 2 x nf matrix
	ans <- lapply(seq.int(nf), function(i)
		affine(x = getFrame(x,i, type = "render"), m = m[,,i], filter = filter, ...))
	if (length(ans) == 1)
		return(ans[[1]])
	else
		return(combine(ans))
}
