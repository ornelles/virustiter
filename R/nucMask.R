#' Generate Nuclear Mask from DNA Image
#'
#' Generate and return an integer Image mask from a fluorescent DNA image
#' such as one obtained by DAPI staining. 
#'
#' @param dna Fluorescent DNA image \emph{or} character vector representing
#'   path(s) to DNA image file(s).
#' @param width Largest nuclear width (diameter) used as width parameter
#'   for \code{thresh2}.
#' @param offset Offset parameter for \code{thresh2}. Use default of 0.05
#'   for typical images, 0.01 for low contrast images.
#' @param size Radius for \code{medianFilter}, integer. Use 2 for routine
#'   images, 0 to skip \code{medianFilter}.
#' @param sigma Standard deviation for \code{gblur}, 2 for routine,
#'   5 for finely detailed.
#' @param radius Radius for code{gblur}, default of 2 * ceiling(3 * sigma)
#'   + 1 for typical imagers, numbers smaller than default of 13 for
#'   smaller features.
#' @param gamma Exponent for \code{DNA^gamma} transformation.
#'
#' @details
#'
#' Generate an integer object mask representing segmented nuclei. The
#' argument can be either a monochrome nuclear fluorescent image of one
#' or more dimensions or a character vector of path(s) to such images.
#'
#' Optimal conditions may need to be found by empirically adjusting the
#' arguments, especially \code{width} and \code{offset}.
#'
#' The image or images will be processed sequentially by (1) gamma
#' transformation, (2) normalization, (3) \code{medianFilter()} with
#' \code{size} (if non-zero), (4) \code{gblur()} with \code{sigma} and 
#' \code{radius}, (5) \code{thresh2()} with \code{width} and \code{offset},
#' (6) \code{fillHull()}, (7) \code{distmap()} and (8) \code{watershed()}.
#'
#' @return
#'
#' A single object holding an integer Image mask for each DNA image.
#'
#' @examples
#' xm <- nucMask(system.file("extdata", "by_folder/b2/file001.tif", package = "virustiter"))
#' max(xm) # total number of nuclei
#' plot(colorLabels(xm))
#'
#' @import EBImage
#'
#' @export
#'
nucMask <- function(dna, width = 36, offset = 0.05, size = 2, sigma = 2,
	radius =  2 * ceiling(3 * sigma) + 1, gamma = 1)
{
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
	x <- gblur(x, sigma, radius)
	x <- thresh2(x, width = width, offset = offset)
	x <- fillHull(x)
	x <- distmap(x)
	return(watershed(x))
}
