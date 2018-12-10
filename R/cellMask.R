#' Voronoi Cell Segmentation from Nuclei
#' 
#' Define boundaries between adjacent cells (regions) from a
#' nuclear mask and optional cellular mask with the
#' \code{\link[EBImage]{propagate}} function. 
#' 
#' @param seeds An \code{Image} object, array, or list of these objects
#'   containing the seeds of identified regions. This would typically
#'   be a segmented nuclear mask or list of masks such as those generated
#'   by \code{nucMask()}.
#' @param mask An optional \code{Image} object, array, or a \code{list}
#'   these objects containing a binary mask defining regions of the image
#'   to be segmented. If this value is \code{NULL}, the nuclear mask
#'   (\code{seeds}) will be expanded for Voronoi segmentation. If the
#'   argument \code{seeds} is a list, \code{mask} must be a similar list
#'   of objects or arrays.
#' @param brush Size of the brush to expand the nuclear mask as an 
#'   odd number of pixels. If this value is \code{NULL}, the mean value of 
#'   the semi-major axis of the nuclei will be used. 
#' @param lambda A numeric value used by \code{propagate()} determining 
#'   the trade-off between the Euclidean distance in the image plane and the 
#'   contribution of the gradient. See \code{\link[EBImage]{propagate}}
#'   for details. 
#'
#' @details
#'
#' A mask to define \emph{approximate} cellular boundaries will be created
#' from a nuclear mask in \code{seeds} and an optional cytoplasmic
#' mask \code{mask}. If the second argument, \code{mask}, is \code{NULL},
#' the nuclear mask will be dilated with a disc-shaped brush of size equal
#' to \code{brush} or, if \code{brush} is \code{NULL}, the average semi-major
#' axis of all nuclei. If \code{mask} is not \code{NULL}, \code{mask} must be
#' a binary mask defining the limits for the Voronoi segmentation based on
#' the seeds provided in \code{seeds}. Such a binary mask can be created from
#' a non-specific cytoplasmic stain such as actin or a diffuse membrane stain.
#'
#' To create a cytoplasmic mask that excludes the nucleus, simply subtract
#' the nuclear mask from the cell mask as shown below. Use \code{erode} or
#' \code{dilate} to adjust the nuclear mask to include more or less of the 
#' peri-nuclear region. 
#'
#' \preformatted{
#'  cmask <- cellMask(nmask) - nmask # when both are single objects
#'  cmask <- cellMask(nmask) - dilate(nmask makeBrush(5, "disc")) # less nucear 
#'  cmask <-lapply(nmask, function(nm) cellMask(nm) - nm) # for list objects
#' }
#'
#' @return
#' 
#' An \code{Image} object produced by \code{propagate()} containing the labeled
#' objects (cells) or a \code{list} of the same.
#'
#' @examples
#'   x <- readImage(system.file("extdata", "by_folder/b2/file001.tif", package = "virustiter"))
#'   y <- readImage(system.file("extdata", "by_folder/b2/file002.tif", package = "virustiter"))
#'   nm <- nucMask(x)
#'   cm <- cellMask(nm)
#'   img <- rgbImage(red = normalize(y) * 0.2, green = normalize(y) * 0.8)
#'   img <- paintObjects(nm, img, col = "cyan")
#'   img <- paintObjects(cm, img, col = "red")
#'   plot(img)
#'
#' @import EBImage
#'
#' @export
#'
cellMask <- function(seeds, mask = NULL, brush = NULL, lambda = 1e-4)
{
# require seeds to be an integer mask or list of the same
	if (is(seeds, "list")) {
		sel <- sapply(seeds, function(x) is.integer(imageData(x)))
		if (!all(sel))
			stop("'", deparse(substitute(seeds)), "' is a list but not all are integer Image masks")
	}
	else if (!is.integer(imageData(seeds)))
		stop("'", deparse(substitute(seeds)), "' is not an integer Image mask")

# process function
	.proc <- function(seeds, mask, brush, lambda)
	{
	# ensure three-dimensions present
		dm <- dim(seeds)
		if (length(dm) == 2)
			dim(seeds) <- c(dm, 1)
	# create mask from seeds if mask == NULL
		if (is.null(mask)) {
			if (is.null(brush))
				brush <- mean(apply(seeds, 3,
					function(x) mean(computeFeatures.moment(x)[,"m.majoraxis"])))
			brush <- 2*round(brush)%/%2 + 1 # make odd
			mask <- dilate(seeds, makeBrush(brush, shape = "disc", step = FALSE))
			mask <- fillHull(mask)
		}
	# ensure that mask is appropriate
		if (!is.integer(imageData(mask)))
			stop("'", deparse(substitute(mask)), "' is not a binary or integer Image mask")
		if (length(dim(mask)) == 2)
			dim(mask) <- c(dim(mask), 1)
		if (!identical(dim(mask), dim(seeds)))
			stop("dimensions of 'seeds' and 'mask' are not same")
	# restore dimensions and return results from propagate
		dim(seeds) <- dim(mask) <- dm
		return(propagate(Image(0, dm), seeds = seeds, mask = mask, lambda = lambda))
	}
# dispatch function accordingly
	if (is(seeds, "Image"))
		ans <- .proc(seeds = seeds, mask = mask, brush = brush, lambda = lambda)
	else
		ans <- lapply(seeds, .proc, mask = mask, brush = brush, lambda = lambda)
	return(ans)
}
