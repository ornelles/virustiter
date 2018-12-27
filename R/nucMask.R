#' Generate Nuclear Mask from DNA Image
#'
#' Generate an integer \code{Image} mask from a fluorescent DNA image. 
#'
#' @param dna Fluorescent DNA \code{Image} or list of fluorescent DNA
#'   \code{Image}s.
#' @param width Maximum nuclear diameter (in pixels) to be used as
#'   the \code{width} parameter for \code{thresh2}.
#' @param offset Offset parameter for \code{thresh2}. Use 0.05
#'   for typical images, use 0.01 for low contrast images.
#' @param size Radius (in pixels) for \code{medianFilter} as an integer.
#'   Use 2 for typical images, use 0 to skip \code{medianFilter}.
#' @param sigma Standard deviation for \code{gblur}, use 2 for typical
#'   images, use 5 for finely detailed images.
#' @param radius Radius for \code{gblur}, use default value of
#'   \code{2 * ceiling(3 * sigma) + 1} for typical images, use numbers
#'   smaller than the default of 13 for images with smaller nuclei.
#' @param gamma Exponent used for \code{DNA^gamma} transformation.
#'
#' @details
#'
#' Generate an integer object mask (or list of masks) representing
#' segmented nuclei. The argument \code{nuc} must be a grayscale nuclear image
#' of one or more dimensions \emph{or} a list of such images.
#'
#' Optimal conditions for detecting and segmented nuclei may require empirically
#' adjusting the arguments, especially \code{width} and \code{offset}.
#'
#' The image or images will be processed sequentially by (1) an optional gamma 
#' transformation, (2) normalization, (3) a \code{medianFilter()} with argument 
#' \code{size} if \code{size} is non-zero, (4) the \code{gblur()} filter with
#' arguments \code{sigma} and \code{radius}, (5) thresholding with \code{thresh2()}
#' with arguments \code{width} and \code{offset}, (6) \code{fillHull()},
#' (7) \code{distmap()}, and (8) \code{watershed()}. 
#'
#' @return
#'
#' An integer \code{Image} mask or list of integer \code{Image} masks.
#'
#' @examples
#'   f.ex <- system.file("extdata", "by_folder/b2/file003.tif",package = "virustiter")
#'   nuc.ex <- readImage(f.ex) # single nuclear image
#'   xm0 <- nucMask(nuc.ex)
#'   plot(colorLabels(xm0))
#'
#' @import EBImage
#'
#' @export
#'
nucMask <- function(dna, width = 36, offset = 0.05, size = 2, sigma = 2,
	radius =  NULL, gamma = 1)
{
	if (missing(dna)) {
		usage <- c("nucMask argument hints:",
			'  dna: fluorescent DNA image or list of images',
			'  width: maximum nuclear diameter',
      '  other argments passed to thresh2(), medianFilter() and gblur()')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}
# internal function to exclude edge objects 
	.edge <- function(v, border) {
		nx <- dim(v)[1]
		ny <- dim(v)[2]
		ans <- unique(c(v[1:border,], v[(nx-border):nx,], v[,1:border], v[,(ny-border):ny]))
		return(ans[ans != 0])
	}

# internal main process function
	.proc <- function(x, width, offset, size, sigma, radius, gamma, border)
	{
	# parameter check and process
		if (gamma != 1) x <- x^gamma
		x <- normalize(x)
		if (!is.null(size) && !is.na(size) && size != 0) {
			size <- as.integer(size)
			x <- medianFilter(x, size)
		}
		if (is.null(radius)) radius <- 2 * ceiling(3 * sigma) + 1
		x <- gblur(x, sigma, radius)
		x <- thresh2(x, width = width, offset = offset)
		x <- fillHull(x)
		x <- distmap(x)
		x <- watershed(x)
		return(x)
	}

# check argument and dispatch function accordingly
	if (is(dna, "Image"))
		ans <- .proc(dna, width = width, offset = offset, size = size,
			sigma = sigma, radius = radius, gamma = gamma)
	else if (all(sapply(dna, is, "Image")))
		ans <- lapply(dna, .proc, width = width, offset = offset, size = size,
			sigma = sigma, radius = radius, gamma = gamma)
	else
		stop("'dna' must be an Image or list of Images")
	return(ans)
}
