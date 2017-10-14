#########################################################################################
# displayPairs
#
# display fluorescent (or DAPI) files in directory with file named in 'f'
#
#########################################################################################

displayPairs <- function(f, dna = FALSE, k.nuc=1.2, width=32, offset=0.05) {
	if(require(EBImage) == FALSE)
		stop("requires EBImage")
	if (missing(f)) f <- file.choose()
	path <- dirname(f)
	ff <- list.files(path, pattern="tif$", full=TRUE, ignore.case = TRUE)
	x <- readImage(ff[seq(1, length(ff), 2)])
	xb <- normalize(x)
	xb <- medianFilter(xb, 2)
	xb <- gblur(xb, 2)
	xt <- thresh(xb, w=width, h=width, offset=offset)
	xt <- fillHull(xt)
	xd <- distmap(xt)
	xw <- watershed(xd)
	y <- readImage(ff[seq(2,length(ff),2)])
	if (dna)
		img <- rgbImage(blue=normalize(x, separate=FALSE))
	else
		img <- rgbImage(green=normalize(y, separate=FALSE))
	img <- paintObjects(xw, img, col="yellow")
	display(img, title=basename(path))
	invisible(img)
}
