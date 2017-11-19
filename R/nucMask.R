#########################################################################################
# nucMask
#
# generate and return nuclear mask from dapi image or dapi file, uses thresh2
#
#	dapi	dapi image or character string of dapi file
#	width	largest nuclear width used as width parameter for thresh2
#	offset	offset parameter for thresh2, default of 0.05, use 0.01 for low contrast
#	size	radius for median filter (integer), 2 for routine images
#	sigma	standard deviation for Gaussian blur, 2 for routine, 5 for finely detailed 
#	gamma	dapi^gamma transformation
#
#########################################################################################

nucMask <- function(dapi, width = 36, offset = 0.05, size = 2, sigma = 2, gamma = 1)
{
	if (is.character(dapi) && file.exists(dapi))
		x <- readImage(dapi)
	else
		x <- dapi
	size <- as.integer(size)
	x <- x^gamma
	x <- normalize(x)
	x <- medianFilter(x, size)
	x <- gblur(x, sigma)
	x <- thresh2(x, width = width, offset = offset)
	x <- fillHull(x)
	x <- distmap(x)
	return(watershed(x))
}
