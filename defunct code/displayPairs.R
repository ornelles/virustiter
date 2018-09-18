#########################################################################################
# displayPairs
#
# display fluorescent (or DAPI) files in directory with file named in 'f'
#
# Arguments
#  f		DAPI file in directory with paired images
#  dnaStain	if TRUE, show masks on blue nuclei otherwise show on 2nd image as green
#  col		overlay color for mask
#  opac		opacity for mask
# Optional arguments passed to nuclearMask
#  width	largest nuclear width used as width parameter for thresh2
#  offset	offset parameter for thresh2, default of 0.05, use 0.01 for low contrast
#  size		radius for median filter (integer), 2 for routine images
#  sigma	standard deviation for Gaussian blur, 2 for routine, 5 for finely detailed 
#  gamma	dapi^gamma transformation
# ...	parameters passed to display
#
# Return
#	side effect of plotting paired images, returned invisibly
#
#########################################################################################

displayPairs <- function(f, dnaStain = FALSE, col = "yellow", opac = 0.5, type = "tiff",
	width = NULL, offset = NULL, size = NULL, sigma = NULL, gamma = NULL, ...)
{
	if(require(EBImage) == FALSE)
		stop("requires EBImage")
	ff <- list.images(dirname(f), type = type)
	maskArgs <- list(width=width, offset=offset, size=size, sigma=sigma, gamma=gamma)
	maskArgs <- maskArgs[!sapply(maskArgs, is.null)]
	x <- readImage(ff[seq(1, length(ff), 2)])
	xw <- do.call(nucMask, c(list(x), maskArgs))
	x <- normalize(x, separate = FALSE)
	y <- readImage(ff[seq(2, length(ff), 2)])
	y <- normalize(y, separate = FALSE)
	if (dnaStain)
		img <- rgbImage(blue = x, green = 0.4 * x)
	else
		img <- rgbImage(green = y, red = 0.4 * y)
	img <- paintObjects(xw, img, col = col, opac = opac)
	display(img, title = basename(f), ...)
	invisible(img)
}
