#' Identify Objects along Edge of Image
#' 
#' Identify objects within \code{border} pixels of the edge of an
#' integer \code{Image} object. 
#' 
#' @param x An \code{Image} object in \code{Grayscale} color mode or an
#'   array containing object masks. Object masks are sets of pixels with
#'   the same unique integer value.
#' @param border Number of pixels from the edge to include.
#' 
#' @return
#'
#' A vector of integers or list of integers identifying objects at the edge
#' 
#' @examples
#'   x <- readImage(system.file("extdata", "by_folder/b4/file001.tif", package = "virustiter"))
#'   nm0 <- nucMask(x)
#'   sel <- edgeObjects(nm0, border = 16)
#'   nm1 <- rmObjects(nm0, sel)
#'   plot(colorLabels(combine(nm0, nm1)), all = TRUE, nx = 1)
#'
#' @import EBImage
#'
#' @export
#'
edgeObjects <- function (x, border = 1)
{
# parameter check
	if (colorMode(x) != 0 || !is.integer(imageData(x)))
		stop("grayscale integer mask needed")
	border <- as.integer(border)

# prepare answer in case that no objects are found
	dm <- dim(x)
	if (length(dm) == 2)
		ans <- numeric(0)
	else
		ans <- rep(list(numeric(0)), dm[3])

# deal with silly border values
	if (border <= 0)
		return(ans)
	if (2*border + 1 > dm[1] || 2*border + 1 > dm[2]) {
		warning("border too large for image size")
		return(ans)
	}

# process function
	.edge <- function(v, border) {
		nx <- dim(v)[1]
		ny <- dim(v)[2]
		ans <- unique(c(v[1:border,], v[(nx-border):nx,], v[,1:border], v[,(ny-border):ny]))
		return(ans[ans!=0])
	}

# dispatch on the dimensions of the argument
	if (length(dm) == 2)
		ans <- .edge(x, border = border)
	else {
		ans <- apply(x, 3, .edge, border = border)
		if (is(ans, "matrix")) # ensure that a list is returned
			ans <- split(ans, c(col(ans)))
	}
	return(ans)
}