#' Generate Nuclear Mask from DNA Image
#'
#' Generate and return an integer \code{Image} mask from a fluorescent DNA image. 
#'
#' @param dna Fluorescent DNA \code{Image} object \emph{or} a character
#'   vector representing a path or paths to DNA image file(s).
#' @param width Largest nuclear width (diameter) used as width parameter
#'   for \code{thresh2}.
#' @param offset Offset parameter for \code{thresh2}. Use 0.05
#'   for typical images, use 0.01 for low contrast images.
#' @param size Radius for \code{medianFilter}, integer. Use 2 for typical
#'   images, use 0 to skip \code{medianFilter}.
#' @param sigma Standard deviation for \code{gblur}, use 2 for typical
#'   images, use 5 for finely detailed images.
#' @param radius Radius for \code{gblur}, use default of 2 * ceiling(3 * sigma)
#'   + 1 for typical images, use numbers smaller than the default of 13 for
#'   images with smaller nuclei.
#' @param gamma Exponent for \code{DNA^gamma} transformation.
#' @param border Exclude objects within this many pixels from the edge.  
#'
#' @details
#'
#' Generate an integer object mask representing segmented nuclei. The
#' argument can be either a monochrome nuclear fluorescent image of one
#' or more dimensions \emph{or} a character vector of path(s) to such images.
#'
#' Optimal conditions may need to be found by empirically adjusting the
#' arguments, especially \code{width} and \code{offset}.
#'
#' The image or images will be processed sequentially by (1) an optional gamma 
#' transformation, (2) normalization, (3) \code{medianFilter()} with argument 
#' \code{size} if \code{size} is non-zero, (4) \code{gblur()} with arguments 
#' \code{sigma} and \code{radius}, (5) \code{thresh2()} with arguments 
#' \code{width} and \code{offset}, (6) \code{fillHull()}, (7) \code{distmap()}, 
#' (8) \code{watershed()} and then (9) optionally removing objects along
#' the edge of the image.
#'
#' @return
#'
#' A single object holding an integer \code{Image} mask for each DNA image.
#' A mask for a single image of dimension \code{dim(dna) = c(nx, ny)} will be
#' returned with a third dimension equal to 1. \code{dim(mask) = c(x, y, 1)}. 
#'
#' @examples
#'   f.example <- system.file("extdata", "by_folder/b2/file003.tif",package = "virustiter")
#'   xm0 <- nucMask(f.example)
#'   max(xm0) # total number of nuclei
#'   xm4 <- nucMask(f.example, border = 4)
#'   max(xm4) # total number of nuclei
#'   opar <- par(mfrow = c(1, 2))
#'   plot(colorLabels(xm0))
#'   plot(colorLabels(xm4))
#'   par(opar)
#'
#' @import EBImage
#'
#' @export
#'
nucMask <- function(dna, width = 36, offset = 0.05, size = 2, sigma = 2,
	radius =  NULL, gamma = 1, border = 0)
{
# internal function to exclude edge objects 
	.edge <- function(v, border) {
		nx <- dim(v)[1]
		ny <- dim(v)[2]
		ans <- unique(c(v[1:border,], v[(nx-border):nx,], v[,1:border], v[,(ny-border):ny]))
		return(ans[ans != 0])
	}
# parameter check and process
	if (is.character(dna) && file.exists(dna))
		x <- readImage(dna)
	else
		x <- dna
	if (gamma != 1)
		x <- x^gamma
	x <- normalize(x)
	if (!is.null(size) && !is.na(size) && size != 0) {
		size <- as.integer(size)
		x <- medianFilter(x, size)
	}
	if (is.null(radius))
		radius <- 2 * ceiling(3 * sigma) + 1
	x <- gblur(x, sigma, radius)
	x <- thresh2(x, width = width, offset = offset)
	x <- fillHull(x)
	x <- distmap(x)
	x <- watershed(x)

	if (border > 0) {
		if (2*border + 1 > dim(x)[1] || 2*border + 1 > dim(x)[2])
			warning("'border' too large for image size")
		else {
			len <- length(dim(x))
			if (len == 2)
				sel <- .edge(x, border)
			else
				sel <- apply(x, len, .edge, border = border)
			x <- rmObjects(x, sel)
		}
	}
	if (length(dim(x)) == 2)
		dim(x) <- c(dim(x), 1)
	return(x)
}
