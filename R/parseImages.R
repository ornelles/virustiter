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
#	display	display colorMap of nuclei via browser
#	nx		nrows to be handed to tile()
#	mask	if TRUE, the return value is a list of image data and nuclear masks
# Arguments passed to nucMask
#	width	largest nuclear width used as w and h parameters for thresh2
#	offset	offset parameter for thresh2, default of 0.05, use 0.01 for low contrast
#	sigma	radius for medianFilter and gblur, default value of 2 for routine images
#
# Returns 	either image image data including well, row, and column and 
#			or a two-element list with value = data and mask = nuclear masks
#
# Use mergePdata on value before using further
#
#########################################################################################

parseImages <- function(dnaFile, k.upper = 3, k.lower = 1.2, width = 36,
	offset = 0.05, sigma = 2, display = TRUE, nx = NULL, mask = FALSE)
{
	library(EBImage)

# examine dnaFile and confirm monochrome
	if (missing(dnaFile)) dnaFile <- file.choose()
	img <- suppressWarnings(readImage(dnaFile))
	if (colorMode(img) != 0) stop("only grayscale images can be used")
	mask <- as.logical(mask)

# determine if image files are a collection of folders or an image stack
# based on dimensions of image in dnaFile
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
				display = display, nx = nx, mask = mask)
	}
	if (mask == TRUE) {
		value <- do.call(rbind, lapply(ret, "[[", 1))
		value <- cbind(dir = factor(dname), value)
		rownames(value) <- NULL
		nmask <- do.call(combine, lapply(ret, "[[", 2))
		ans <- list(value = value, mask = nmask)
	}
	else {
		value <- do.call(rbind, ret)
		value <- cbind(dir = factor(dname), value)
		rownames(value) <- NULL
		ans <- value
	}
	if (!display) close(pb)
	return(ans)
}
#
# process paired DAPI and fluorescent images in named files
#
.processByFolder <- function(imageFiles, k.upper, k.lower, width,
		offset, sigma, display, nx, mask)
{
	img <- suppressWarnings(readImage(imageFiles))
	dapi <- img[,,seq(1, dim(img)[3], 2)]	# must have 3 dimensions
	flr <- img[,,seq(2, dim(img)[3], 2)]
	if (is.null(nx))
		nx <- ceiling(sqrt(dim(dapi)[3]))
	xw <- nucMask(dapi, sigma = sigma, width = width, offset = offset)

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
	if (mask == TRUE)
		return(list(value = res, mask = xw))
	else
		return(res)
}

#
# process alternately stacked DAPI and fluorescent images in named file
#
.processByStack <- function(imageStack, k.upper, k.lower, width,
		offset, sigma, display, nx, mask)
{
	img <- suppressWarnings(readImage(imageStack))
	dapi <- img[,,seq(1, dim(img)[3], 2)]	# must have 3 dimensions
	flr <- img[,,seq(2, dim(img)[3], 2)]
	if (is.null(nx))
		nx <- ceiling(sqrt(dim(dapi)[3]))
	xw <- nucMask(dapi, sigma = sigma, width = width, offset = offset)

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
	if (mask == TRUE)
		return(list(value = res, mask = xw))
	else
		return(res)
}
