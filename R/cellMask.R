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
#'   first argument \code{seeds} is a list, \code{mask} must be a similar
#'   list of objects or arrays.
#' @param brush Size of the brush to expand the nuclear mask as an 
#'   odd number of pixels. If this value is \code{NULL}, the mean value of 
#'   the semi-major axis of the nuclei will be used. If this value is
#'   negative, the resulting mask will be dilated by the brush size. 
#'   IN PROGRESS - NOT YET.
#' @param lambda A numeric value used by \code{propagate()} determining 
#'   the trade-off between the Euclidean distance in the image plane and the 
#'   contribution of the gradient. See \code{\link[EBImage]{propagate}}
#'   for details. 
#'
#' @details
#'
#' A mask to define \emph{approximate} cellular boundaries will be created
#' from a nuclear mask in \code{seeds} and an optional cytoplasmic
#' mask \code{mask}. If the options second argument \code{mask} is \code{NULL},
#' the nuclear mask will be dilated with a disc-shaped brush of size equal
#' to \code{brush} or, if \code{brush} is \code{NULL}, the average semi-major
#' axis of all nuclei. If \code{mask} is not \code{NULL}, \code{mask} must be
#' a binary (or integer) mask defining the limits for the Voronoi segmentation
#' based on the seeds provided in \code{seeds}. Such a binary mask can be
#' created by thresholding a non-specific widespread cytoplasmic signal such
#' antibody labeling for actin or a diffuse membrane stain.
#'
#' To create a cytoplasmic mask that excludes the nucleus, simply subtract
#' the nuclear mask from the cell mask as shown below. Use \code{erode()} or
#' \code{dilate()} to adjust the nuclear mask to include more or less of the 
#' peri-nuclear region.
#'
#' To create a smaller nuclear mask, use a negative brush value.
#'
#' \preformatted{
#' # When both are single objects
#'   cytoplasm <- cellMask(nmask) - nmask
#'   cytoplasm2 <- cellMask(nmask) - cellMask(nmask, 5) # capture less of nucleus 
#'   small.nmask <- cellMask(nmask, brush = -5)
#'
#' # When both are list objects
#'   cmask <-lapply(nmask, function(nm) cellMask(nm) - nm) # for list objects
#' }
#'
#' @return
#' 
#' An \code{Image} object produced by \code{propagate()} containing the labeled
#' objects (cells) or a \code{list} of the same.
#'
#' @examples
#'   x <- readImage(system.file("extdata", "by_folder/a4/file001.tif", package = "virustiter"))
#'   y <- readImage(system.file("extdata", "by_folder/a4/file002.tif", package = "virustiter"))
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

# if mask is used, require it to be an integer mask or list of the same
	if (!is.null(mask)) {
		if (is(mask, "list")) {
			sel <- sapply(mask, function(x) is.integer(imageData(x)))
			if (!all(sel))
				stop("'", deparse(substitute(mask)), "' is a list but not all are integer Image masks")
		}
		else if (!is.integer(imageData(mask)))
			stop("'", deparse(substitute(mask)), "' is not an integer Image mask")
	}

# if brush is present, ensure that it is an integer
	if (!is.null(brush))
		brush <- as.integer(brush)

# process function
	.proc <- function(seeds, mask, brush, lambda)
	{
	# ensure three-dimensions present
		dm <- dim(seeds)
		if (length(dm) == 2)
			dim(seeds) <- c(dm, 1)
	# create mask from seeds if mask == NULL
		if (is.null(mask)) {
			if (is.null(brush)) {
				brush <- mean(apply(seeds, 3,
					function(x) mean(computeFeatures.moment(x)[,"m.majoraxis"])))
				brush <- as.integer(brush)
			}
		# apply erosion or dilation
			brush <- 2*(brush + ifelse(brush < 0, -1, 0))%/%2 + 1 # ensure odd number
			if (brush < 0)
				mask <- erode(seeds, makeBrush(-brush, "disc"))
			else if (brush > 0)
				mask <- dilate(seeds, makeBrush(brush, "disc"))
			mask <- fillHull(mask)
			dim(mask) <- dm
		}
	# ensure that mask is appropriate
		if (!is.integer(imageData(mask)))
			stop("'", deparse(substitute(mask)), "' is not a binary or integer Image mask")
		if (!identical(dim(mask), dm))
			stop("dimensions of 'seeds' and 'mask' are not same")
	# restore dimensions and return results from propagate
		dim(seeds) <- dm
		return(propagate(Image(0, dm), seeds = seeds, mask = mask, lambda = lambda))
	}

# dispatch function accordingly
	if (is(seeds, "Image"))
		ans <- .proc(seeds = seeds, mask = mask, brush = brush, lambda = lambda)
	else if (is(seeds, "list") && is.null(mask))
		ans <- lapply(seeds, function(s) .proc(s, mask = NULL, brush = brush, lambda = lambda))
	else if (is(seeds, "list") && is(mask, "list") && length(seeds) == length(mask))
		ans <- Map(function(s, m) .proc(s, m, brush = brush, lambda = lambda),
			seeds, mask)
	else
		stop(deparse(substitute(seeds)), " and ", deparse(substitute(mask)),
			" are not compatible")
	return(ans)
}
