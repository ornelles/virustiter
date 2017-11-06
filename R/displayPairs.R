#########################################################################################
# displayPairs
#
# display fluorescent (or DAPI) files in directory with file named in 'f'
#
# Arguments
# f			path to DAPI image in directory with paired images
# dnaStain	if TRUE, add blue nuclei to image, else show 2nd as green
# width, offset, sigma	paramaters for nucMask
# col		overlay color for mask 
# ...	parameters passed to display
#
# Return
#	side effect of plotting paired images, returned invisibly
#
#########################################################################################

displayPairs <- function(f, dnaStain = FALSE, width = 32,
	offset = 0.05, sigma = 2, col = "lightyellow", ...)
{
	if(require(EBImage) == FALSE)
		stop("requires EBImage")
	if (missing(f)) f <- file.choose()
	path <- dirname(f)
	ff <- list.files(path, pattern="tif$", full=TRUE, ignore.case = TRUE)
	x <- readImage(ff[seq(1, length(ff), 2)])
	xw <- nucMask(x, width=width, offset=offset, sigma=sigma)
	y <- readImage(ff[seq(2, length(ff), 2)])
	if (dnaStain)
		img <- rgbImage(blue = normalize(x, separate = FALSE))
	else
		img <- rgbImage(green = normalize(y, separate = FALSE))
	img <- paintObjects(xw, img, col = col)
	display(img, title = basename(path), ...)
	invisible(img)
}
