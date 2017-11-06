#########################################################################################
# nucMask
#
# generate and return nuclear mask from dapi image or dapi file, uses thresh2
#
#	dapi	dapi image or character string of dapi file
#	width	largest nuclear width used as w and h parameters for .thresh2()
#	offset	offset parameter for .thresh2(), default of 0.05, use 0.01 for low contrast
#	sigma	radius for medianFilter and gblur, default value of 2 routine 
#
#########################################################################################

nucMask <- function(dapi, width = 36, offset = 0.05, sigma = 2) {
	if (is.character(dapi) && file.exists(dapi))
		dapi <- readImage(dapi)
	x <- normalize(dapi)
	x <- medianFilter(x, sigma)
	x <- gblur(x, sigma)
	x <- thresh2(x, w = width, offset = offset)
	x <- fillHull(x)
	x <- distmap(x)
	return(watershed(x))
}

#
# adaptive threshold using disc instead of box
#
thresh2 <- function(x, w, offset) {
	r <- w - w%%2 + 1
	f <- makeBrush(r, shape="disc")
	f <- f/sum(f)
	return (x > (filter2(x, f) + offset))
}
