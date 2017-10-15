#########################################################################################
# parseImages
#
# Use EBImage tools to analyze image files and replace what the ImageJ code,
# "Multiwell Fluorescent Cell Count," does. Creates similar output to readIJResults().
# Merge this output with a phenotype data file to create a file of the same form.
#
# Arguments
#	dnaFile	optional path to DAPI file under nested folders such as a3/file001.tif
#	k.upper	multiplier for mad value cutoff for nuclear area: mean + k.upper*mad(area)
#	k.lower	multiplier for mad value cutoff for nuclear area: mean - k.lower*mad(area)
#	width	w and h parameters for thresh()
#	offset	offset parameter for thresh(), default of 0.05, use 0.01 for low contrast
#	sigma	radius for medianFilter and gblur use 2 for routine, 
#	draw	draw colorMap of nuclei
#	nx		nrows to be handed to tile()
#
# Returns invisible raw data with extracted well, row, column use mergeData to assign
# additional information before handing to other functions
#
#########################################################################################

parseImages <- function(dnaFile, k.upper = 3, k.lower = 1.2, width = 32, offset = 0.05,
		sigma = 2, draw = TRUE, nx = NULL)
{
	library(EBImage)

# working function to apply to paired dapi and fluorescent images
	.fun <- function(ff, k.upper, k.lower, width, offset, sigma, draw, nx) {
		img <- readImage(ff)
		dapi <- img[,,seq(1, dim(img)[3], 2)]
		flr <- img[,,seq(2, dim(img)[3], 2)]
		if (is.null(nx))
			nx <- ceiling(sqrt(dim(dapi)[3]))
		xb <- normalize(dapi)
		xb <- medianFilter(xb, sigma)
		xb <- gblur(xb, sigma)
		xb <- thresh2(xb, w = width, offset = offset)
		xb <- fillHull(xb)
		xd <- distmap(xb)
		xw <- watershed(xd)
		well <-  well.info(basename(dirname(ff)))$well[1]
		column <- well.info(well)$column
		row <- well.info(well)$row

		fname <- paste(well, basename(ff)[seq(2, length(ff), 2)], sep="/")

		nframes <- dim(xw)[3]
		area <- lapply(1:nframes, function(i) computeFeatures.shape(xw[,,i])[,1])
		xbar <- mean(unlist(area))
		xmad <- mad(unlist(area))
		small <- lapply(area, function(z) which(z < xbar - k.lower*xmad))
		large <- lapply(area, function(z) which(z > xbar + k.upper*xmad))
		xw <- rmObjects(xw, small, reenumerate = FALSE)
		xw <- rmObjects(xw, large)
		if (draw)
			display(tile(colorLabels(xw), nx = nx, lwd = 3, fg.col="white"), title = well)
		res <- data.frame()
		for (i in seq_len(nframes)) {
			area <- computeFeatures.shape(xw[,,i])[,1]
			XY <- computeFeatures.moment(xw[,,i])[,1:2]
			dna <- computeFeatures.basic(xw[,,i], dapi[,,i])[,1]
			val <- computeFeatures.basic(xw[,,i], flr[,,i])[,1]
			res <- rbind(res, data.frame(well = well, column = column, row = row,
				file = fname[i], xm = XY[,1], ym = XY[,2], area, dna, val))
		}
		rownames(res) <- NULL
		return(res)
	}

# extract list of paths to image pairs from dnaFile
	if (missing(dnaFile))
		dnaFile <- file.choose()
	path <- dirname(dirname(dnaFile))
	dname <- basename(path)
	ff <- list.files(path, full = TRUE, recursive = TRUE, pattern = "tif$",
			ignore.case = TRUE)
	g <- basename(dirname(ff))
	fl <- split(ff, g)
	bad <- which(lengths(fl)%%2 != 0)	# check for mismatch images
	if (length(bad))
		stop("odd number of images in: ", paste(names(bad)))

# apply the working function to extract information
	ret <- rep(list(NULL), length(fl))
	if (!draw) pb <- txtProgressBar(min = 1, max = length(fl), style = 3)
	for (i in seq_along(fl)) {
		if (!draw) setTxtProgressBar(pb, i)
		ret[[i]] <- .fun(fl[[i]], k.upper = k.upper, k.lower = k.lower,
				width = width, offset = offset, sigma = sigma, draw = draw, nx = nx)
	}
	ret <- do.call(rbind, ret)
	ret$dname <- factor(dname)
	rownames(ret) <- NULL
	if (!draw) close(pb)
	return(ret)
}
#
# adaptive threshold using disc instead of box
#
thresh2 <- function(x, w=5, offset=0.01) {
	r <- w - w%%2 + 1
	f <- makeBrush(r, shape="disc")
	f <- f/sum(f)
	return (x > (filter2(x, f) + offset))
}
