#' Normalize Image Background
#' 
#' Adjust values of an image to a common background.
#' 
#' @param img Grayscale \code{Image} object.
#' @param inputRange A numeric vector of 2 values specifying the expected
#'   range of intensity values.
#' @param quant Quantile to serve as the common baseline, default of 0.025.
#' @param nonzero Should zero values be ignored in quantile calculation?
#' @param ... Additional arguments are accepted but ignored in order
#'   to remain compatible with previous versions of this function.
#' 
#' @details
#' 
#' Each frame of an \code{Image} will be linearly scaled to have an
#' identical value at the specified quantile. If \code{inputRange} is missing,
#' it will be taken as the actual range of values in \code{img}.
#' Each frame of the image will be adjusted by subtracting the median value at
#' the specified quantile and adding \code{quant * inputRange[2]}. Negative 
#' values will be clipped to 0. \emph{Note that if \code{inputRange}
#' is not specified, the different image objects will not necessarily share
#' the same background value.} 
#' 
#' @return
#' 
#' \code{Image} object linearly scaled to have a common background.
#' 
#' @examples
#'   cell <- readImage(system.file("extdata", "by_folder/b2/file002.tif", 
#'    package = "virustiter"))
#'   N <- prod(dim(cell)[1:2])
#'   cells <- Image(rep(imageData(cell), 4), dim = c(dim(cell), 4))
#'   cells[,,2] <- cell + 0.005
#'   cells[,,3] <- cell + 0.01
#'   cells[,,4] <- cell + 0.02
#'   cells <- normalize(cells, separate = FALSE)
#'   plot(cells, all = TRUE)
#'   plot(bnormalize(cells), all = TRUE)
#' 
#' @import EBImage
#' 
#' @export
#' 
bnormalize <- function(img, inputRange, quant = 0.025, nonzero = FALSE, ...)
{
	if (!is(img, "Image") || colorMode(img) != 0)
		stop("grayscale image required")
	if (quant < 0 | quant > 1)
		stop("'quant' must be between 0 and 1")
	if (missing(inputRange))
		inputRange <- range(img)
	base <- quant * inputRange[2]

	dm <- dim(img)
	if (length(dm) == 2) {
		if (nonzero == TRUE)
			setpoint <- quantile(img[img > 0], quant)
		else
			setpoint <- quantile(img, quant)
		img <- img - setpoint + base
	}
	else {
		if (nonzero == TRUE)
			setpoint <- apply(img, 3, function(v) quantile(v[v > 0], quant))
		else
			setpoint <- apply(img, 3, quantile, quant)
		bgndImg <- Image(rep(setpoint, each = prod(dm[1:2])), dim = dm)
		img <- img - bgndImg + base
	}
	img[img < 0] <- 0
	return(img)
}
