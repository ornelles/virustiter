#########################################################################################
# nucMask
#
# generate and return nuclear mask from dapi image or dapi file, uses thresh2
#
#	dapi	dapi image or character string of dapi file
#	width	largest nuclear width used as width parameter for thresh2
#	offset	offset parameter for thresh2, default of 0.05, use 0.01 for low contrast
#	size	integer, radius for median filter, 2 for routine images,
#			0 to exclude medianFilter
#	sigma	standard deviation for Gaussian blur, 2 for routine, 5 for finely detailed
#	radius	radius for gblur, default of 2 * ceiling(3 * sigma) + 1
#			use smaller numbers than default of 13 for smaller features
#	gamma	dapi^gamma transformation
#
# update 6/2018
#	add 'radius' to pass to gblur() use small values for small features
#	accept values of 0 for 'size' to exclude medianFilter()
#########################################################################################

nucMask <- function(dapi, width = 36, offset = 0.05, size = 2, sigma = 2,
	radius =  2 * ceiling(3 * sigma) + 1, gamma = 1)
{
	if (is.character(dapi) && file.exists(dapi))
		x <- readImage(dapi)
	else
		x <- dapi
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
