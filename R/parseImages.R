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
#	display	display colorMap of nuclei
#	nx		nrows to be handed to tile()
#
# Returns raw data with extracted well, row, column. Use mergePdata before using further
#
#########################################################################################

parseImages <- function(dnaFile, k.upper = 3, k.lower = 1.2, width = 32,
	offset = 0.05, sigma = 2, display = TRUE)
{
	library(EBImage)

# examine dnaFile and confirm monochrome
	if (missing(dnaFile)) dnaFile <- file.choose()
	img <- suppressWarnings(readImage(dnaFile))
	if (colorMode(img) != 0) stop("only grayscale images can be used")

# process based on dimensions of image in dnaFile
	if (length(dim(img)) == 2 || dim(img)[3] == 1) { # folders
		path <- dirname(dirname(dnaFile))
		dname <- basename(path)
		ff <- list.files(path, full = TRUE, recursive = TRUE,
				pattern = "tif{1,2}$", ignore.case = TRUE)
		ffsplit <- split(ff, basename(dirname(ff)))
		nff <- length(ffsplit)
		bad <- which(lengths(ffsplit)%%2 != 0)	# mismatched files?
		if (length(bad))
			stop("odd number of images in: ", paste(names(bad)))
		FUN <- .processByFolder
	}
	else {  # image stacks
		path <- dirname(dnaFile)
		dname <- basename(path)
		ff <- list.files(path, full = TRUE, recursive = TRUE,
				pattern = "tif$", ignore.case = TRUE)
		ffsplit <- split(ff, basename(ff))
		nff <- length(ff)
		FUN <- .processByStack
	}
# apply the working function to extract information
	ret <- rep(list(NULL), nff)
	if (!display)
		pb <- txtProgressBar(min = 1, max = length(ffsplit), style = 3)
	for (i in seq_along(ffsplit)) {
		if (!display) setTxtProgressBar(pb, i)
		ret[[i]] <- FUN(ffsplit[[i]], k.upper = k.upper, k.lower = k.lower,
				width = width, offset = offset, sigma = sigma,
				display = display)
	}
	ret <- do.call(rbind, ret)
	ret <- cbind(directory = factor(dname), ret)
	rownames(ret) <- NULL
	if (!display) close(pb)
	return(ret)
}
#
# adaptive threshold using disc instead of box
#
.thresh2 <- function(x, w=5, offset=0.01) {
	r <- w - w%%2 + 1
	f <- makeBrush(r, shape="disc")
	f <- f/sum(f)
	return (x > (filter2(x, f) + offset))
}
#
# process paired DAPI and fluorescent images in named files
#
.processByFolder <- function(imageFiles, k.upper, k.lower, width,
		offset, sigma, display)
{
	img <- suppressWarnings(readImage(imageFiles))
	dapi <- img[,,seq(1, dim(img)[3], 2)]	# must have 3 dimensions
	flr <- img[,,seq(2, dim(img)[3], 2)]
	nx <- ceiling(sqrt(dim(dapi)[3]))
	xb <- normalize(dapi)
	xb <- medianFilter(xb, sigma)
	xb <- gblur(xb, sigma)
	xb <- .thresh2(xb, w = width, offset = offset)
	xb <- fillHull(xb)
	xd <- distmap(xb)
	xw <- watershed(xd)

# add well descriptors
	well <-  well.info(basename(dirname(imageFiles)))$well[1]
	column <- well.info(well)$column
	row <- well.info(well)$row

	fname <- paste(well, basename(imageFiles)[seq(2, length(imageFiles), 2)], sep = "/")

	nframes <- dim(xw)[3]
	area <- lapply(1:nframes, function(i) computeFeatures.shape(xw[,,i])[,1])
	xbar <- mean(unlist(area))
	xmad <- mad(unlist(area))
	small <- lapply(area, function(z) which(z < xbar - k.lower*xmad))
	large <- lapply(area, function(z) which(z > xbar + k.upper*xmad))
	xw <- rmObjects(xw, small, reenumerate = FALSE)
	xw <- rmObjects(xw, large)
	if (display)
		display(tile(colorLabels(xw), nx = nx, lwd = 3, fg.col="white"), title = well)
	res <- data.frame()
	for (i in seq_len(nframes)) {
		area <- computeFeatures.shape(xw[,,i])[,1]
		XY <- computeFeatures.moment(xw[,,i])[,1:2]
		dna <- computeFeatures.basic(xw[,,i], dapi[,,i])[,1]
		val <- computeFeatures.basic(xw[,,i], flr[,,i])[,1]
		res <- rbind(res, data.frame(file = fname[i], well = well, column = column,
			row = row, xm = XY[,1], ym = XY[,2], area, dna, val))
	}
	rownames(res) <- NULL
	return(res)
}

#
# process alternately stacked DAPI and fluorescent images in named file
#
.processByStack <- function(imageStack, k.upper, k.lower, width,
		offset, sigma, display)
{
	img <- suppressWarnings(readImage(imageStack))
	dapi <- img[,,seq(1, dim(img)[3], 2)]	# must have 3 dimensions
	flr <- img[,,seq(2, dim(img)[3], 2)]
	nx <- ceiling(sqrt(dim(dapi)[3]))
	xb <- normalize(dapi)
	xb <- medianFilter(xb, sigma)
	xb <- gblur(xb, sigma)
	xb <- .thresh2(xb, w = width, offset = offset)
	xb <- fillHull(xb)
	xd <- distmap(xb)
	xw <- watershed(xd)

# extract file name and number of frames
	fname <- basename(imageStack)
	nframes <- dim(xw)[3]

	area <- lapply(1:nframes, function(i) computeFeatures.shape(xw[,,i])[,1])
	xbar <- mean(unlist(area))
	xmad <- mad(unlist(area))
	small <- lapply(area, function(z) which(z < xbar - k.lower*xmad))
	large <- lapply(area, function(z) which(z > xbar + k.upper*xmad))
	xw <- rmObjects(xw, small, reenumerate = FALSE)
	xw <- rmObjects(xw, large)
	if (display)
		display(tile(colorLabels(xw), nx = nx, lwd = 3, fg.col="white"), title = fname)
	res <- data.frame()
	for (i in seq_len(nframes)) {
		area <- computeFeatures.shape(xw[,,i])[,1]
		XY <- computeFeatures.moment(xw[,,i])[,1:2]
		dna <- computeFeatures.basic(xw[,,i], dapi[,,i])[,1]
		val <- computeFeatures.basic(xw[,,i], flr[,,i])[,1]
		res <- rbind(res, data.frame(file = fname, frame = i, 
			xm = XY[,1], ym = XY[,2], area, dna, val))
	}
	rownames(res) <- NULL
	return(res)
}
