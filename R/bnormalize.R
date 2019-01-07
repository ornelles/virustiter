#' Normalize Image Background
#' 
#' Adjust images to have a common background value by simple linear shift.
#' 
#' @param img \code{Image} grayscale image object
#' @param size Radius for \code{medianFilter}, integer. Default value of 2 
#'   for typical images; use 0 to skip median filter.
#' @param sigma Standard deviation for \code{gblur}. Default value of 2 
#'   for typical images, use 5 for finely detailed images; use 0 to skip 
#'   Gaussian blur.
#' @param quant Quantile to evaluate as baseline for background signal
#'  (default value of 0.05, use 0.5 for median value).
#' @param offset Value (default of 0.05) added to image after subtracting 
#'  baseline value. 
#' 
#' @details
#' 
#' This is meant to be applied to images with a greater number of pixels 
#' representing background values than foreground values. The baseline
#' for each frame will be determined as \code{quantile(x, quant)} where
#' \code{x} is the pixel intensity. Each image will be adjusted 
#' to have the identical background value at this quantile. The image(s) will
#' first be smoothed by \code{medianFilter()} with radius \code{size} then
#' blurred by \code{gblur} with \code{sigma}. Each processed image will be 
#' adjusted by subtracting the background value and adding \code{offset}. Values
#' less than 0 will be changed to 0.
#' 
#' @return
#' 
#' \code{Image} object of the same size that has been smoothed, blurred 
#' and linearly adjusted.
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
bnormalize <- function(img, size = 2, sigma = 2, quant = 0.05, offset = 0.05)
{
	if (!is(img, "Image") || colorMode(img) != 0)
		stop("grayscale image required")
	if (size > 0) {
		size <- as.integer(size)
		img <- medianFilter(img, size)
	}
	if (sigma > 0)
		img <- gblur(img, sigma)
	if (quant[1] > 1 | quant[1] < 0)
		stop("'quant' must be between 0 and 1")

	dm <- dim(img)
	if (length(dm) == 2) {
		bgnd <- quantile(img, quant)
		img <- img - bgnd + offset
	}
	else {
		bgnd <- apply(img, 3, quantile, quant)
		bgndImg <- Image(rep(bgnd, each = prod(dm[1:2])), dim = dm)
		img <- img - bgndImg + offset
	}
	img[img < 0] <- 0
	return(img)
}
