#########################################################################################
# cellMask
#
# generate and return a cell mask from cell and nuclear seed images using Otsu's
# method for bimodal population on a low-pass filtered cell image
#
#	cell	cell (cytoplasm) image or character string of dapi file
# seeds	segmented nuclear mask, from nucMask for example
#	lpf.width	diameter for low-pass filter
#	thresh	optional explicit threshold if Otsu's method not used
#	brush	brush size for lpf and opening
#	lambda	parameter for propagate
#	bitDepth	bit depth of image for Otsu calculations
#
#########################################################################################

cellMask <- function(cell, seeds, lpf.width = 21, thresh = NULL, brush = 5,
	lambda = 1e-4, bitDepth = 8)
{
	if (is.character(cell) && file.exists(cell))
		x <- readImage(cell)
	else
		x <- cell
	if (!identical(dim(x), dim(seeds)))
		stop("dimensions of 'cell' and 'seeds' are not same")
	cb <- normalize(x)
	w <- 2 * (lpf.width %/% 2) + 1
	filt <- makeBrush(w, shape = "disc", step = FALSE)
	filt <- filt/sum(filt)
	cb <- filter2(cb, filt)
	if (!is.null(thresh))
		ct <- opening(cb > thresh, makeBrush(brush, shape = "disc"))
	else {
		ct <- array(NA, dim = dim(cb))
		thresh <- otsu(cb, levels = 2^bitDepth)
		for (i in seq_along(thresh))
			ct[,,i] <- cb[,,i] > thresh[i]
		ct <- dilate(ct, makeBrush(lpf.width, shape = "Gaussian"))
		ct <- fillHull(ct)
	}
	cmask <- propagate(x, seeds = seeds, mask = ct, lambda = lambda)
}
