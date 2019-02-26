#' Interactively Measure Point-to-Point Distances
#' 
#' Use the \code{locator} function in base graphics to interactively
#' measure distances between pairs of points on an existing image.
#' End the interaction pressing the \strong{Esc} key \emph{or}
#' control-clicking.
#' 
#' @param n	Maximum number of pairs to measure, integer.
#' @param col,type,pch Values handed to \code{locator()}.
#' @param ... Additional values handed to \code{locator()}.
#' 
#' @return
#' 
#' Distances between points with locator in pixels. 
#'
#' @export
#'
p2p <- function(n = 512, col = "magenta", type = "o", pch = 3, ...)
{
	ans <- numeric()
	while (n > 0) {
		p <- locator(2, type = type, pch = pch, col = col, ...)
		if (is.null(p))
			break
		ans <- c(ans, sqrt(sum(sapply(p, diff)^2)))
		n <- n - 1
	}
	return(ans)
}
