#########################################################################################
# trimMask
#
# remove small and large objects from mask
#
#	mask	mask object
#	cutoff	integer, upper and lower absolute cutoff values 
#	k		real, upper and lower multiplier to determine cutoff values
#			if cutoff is NULL from mean (xbar) and MAD (xmad) area
#			cutoff <- c(xbar - k[1] * xmad, xbar + k[2] * xmad)
#			reenumerate	logical, reenumerate before returning 
#
#########################################################################################
trimMask <- function(mask, cutoff = NULL, k = c(1.5, 3), reenumerate = TRUE)
{
	require(EBImage)
	dm <- dim(mask)
	if(length(dm) < 2 || length(dm) > 3)
		stop("'mask' must be a 2- or 3-dimension integer array")
	if (length(dm) == 2)
		dim(mask) <- c(dm, 1)
	nframes <- dim(mask)[3]
	area <- lapply(1:nframes, function(i) computeFeatures.shape(mask[,,i])[,1])
	xbar <- mean(unlist(area))
	xmad <- mad(unlist(area))
	if (is.null(cutoff)) {
		lower <- xbar - k[1] * xmad
		upper <- xbar + k[2] * xmad
	}
	else {
		lower <- max(cutoff[1], min(unlist(area)))
		upper <- min(cutoff[2], max(unlist(area)))
	}
	small <- lapply(area, function(z) which(z < lower))
	large <- lapply(area, function(z) which(z > upper))
	mask <- rmObjects(mask, small, reenumerate = FALSE)
	mask <- rmObjects(mask, large)
	dim(mask) <- dm
	if (reenumerate)
		mask <- reenumerate(mask)
	return(mask)
}
