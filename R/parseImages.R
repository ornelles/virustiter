#########################################################################################
# parseImages
#
# Use EBImage tools to analyze image files and replace what the ImageJ code,
# "Multiwell Fluorescent Cell Count," does. Creates similar output to readIJResults().
# Merge this output with a phenotype data file to create a file of the same form.
#
# Arguments
#	dnaFile	path to DAPI file under nested folders or multilayer tif (A1/file001.tif)
#	display	display color map of nuclear mask with browser
#	mask	if TRUE, the return value is a named list with "data" and "mask"
#			otherwise just the data frame of results are returned
#	pattern	optional grep pattern to select image files
#	type	type of image files as "tiff", "jpeg" or "png"	
# Arguments passed to trimMask
#	cutoff	integer, upper and lower absolute cutoff values 
#	k		real, upper and lower multiplier to determine cutoff values
#			if cutoff is NULL from mean (xbar) and MAD (xmad) area
#			lim <- c(xbar - k[1] * xmad, xbar + k[2] * xmad)
# Arguments passed to nucMask
#	width	largest nuclear width used as width parameter for thresh2
#	offset	offset parameter for thresh2, default of 0.05, use 0.01 for low contrast
#	size	integer, radius for median filter, 2 for routine images,
#			0 to exclude medianFilter
#	sigma	standard deviation for Gaussian blur, 2 for routine, 5 for finely detailed
#	radius	radius for gblur, default of 2 * ceiling(3 * sigma) + 1
#			use smaller numbers than default of 13 for smaller features
#	gamma	exponent for gamma transformation (image^gamma)
#
# Returns 	either parsed image data including well, row, and column and 
#			or a two-element named list with data = parsed imaged data
#			and mask = nuclear masks
#
# Use mergePdata on data before proceeding further with virus titer
#
#########################################################################################

parseImages <- function(dnaFile, cutoff = NULL, k = c(1.2, 3), width = 36,
	offset = 0.05, size = 2, sigma = 2, radius =  2 * ceiling(3 * sigma) + 1,
	gamma = 1, display = TRUE, mask = FALSE, type = "tiff", pattern = NULL)
{
# requires EBImage, ensure appropriate values for parameters 
	library(EBImage)
	mask <- as.logical(mask)
	size <- as.integer(size)

# examine dnaFile and confirm that it is monochrome
	img <- suppressWarnings(readImage(dnaFile))
	if (colorMode(img) != 0) stop("only grayscale images can be used")

# use the dimensions of image in dnaFile to determine if image files are a collection
# of folders or an image stack
	if (length(dim(img)) == 2 || dim(img)[3] == 1) { # folders
		FUN <- .processByFolder
		dname <- basename(dirname(dirname(dnaFile)))
		ff <- list.images(dnaFile, enclosing = 2, type = type, pattern = pattern)
		ffsplit <- split(ff, basename(dirname(ff)))
		nff <- length(ffsplit)
		bad <- which(lengths(ffsplit)%%2 != 0)	# mismatched files?
		if (length(bad))
			stop("odd number of images in: ", paste(names(bad)))
	}
	else {  # image stacks
		FUN <- .processByStack
		dname <- basename(dirname(dnaFile))
		ff <- list.images(dnaFile, enclosing = 1, type = type, pattern = pattern)
		ffsplit <- split(ff, basename(ff))
		nff <- length(ff)
	}

# initialize variable to collect results, determine if progress bar is needed
	ret <- rep(list(NULL), nff)
	if (display == TRUE | length(ffsplit) < 2)
		showProgress <- FALSE
	else {
		pb <- txtProgressBar(min = 1, max = length(ffsplit), style = 3)
		showProgress <- TRUE
	}

# apply the working function to extract information from each set of files
	for (i in seq_along(ffsplit)) {
		if (showProgress)
			setTxtProgressBar(pb, i)
		ret[[i]] <- FUN(ffsplit[[i]], cutoff = cutoff, k = k,
				width = width, offset = offset, size = size, sigma = sigma,
				radius = radius, gamma = gamma)
		if (display) {
			if ("well" %in% names(ret[[i]]$data))
				title <- levels(ret[[i]]$data$well)[1]
			else
				title <- levels(ret[[i]]$data$file)[1]
			nx <- ceiling(sqrt(dim(ret[[i]]$mask)[3]))
			display(tile(colorLabels(ret[[i]]$mask), nx = nx, lwd = 3,
					fg.col="white"), title = title)
		}
	}

# FUN returned the data and mask as a two-element list
# add directory information to data
	value <- do.call(rbind, lapply(ret, "[[", 1))
	value <- cbind(dir = factor(dname), value)
	rownames(value) <- NULL

# if required, convert the list of nuclear masks into a single object
	if (mask == TRUE) {
		nmask <- do.call(combine, lapply(ret, "[[", 2))
		ans <- list(data = value, mask = nmask)
	}
	else
		ans <- value
	if (showProgress)
		close(pb)
	return(ans)
}
#
# process paired DAPI and fluorescent images in named files
# return data.frame as 'data' and nuclear mask as 'mask'
#
.processByFolder <- function(imageFiles, cutoff, k, width,
		offset, size, sigma, radius, gamma)
{
# read image files
	img <- suppressWarnings(readImage(imageFiles))
	dapi <- img[,,seq(1, dim(img)[3], 2)]
	flr <- img[,,seq(2, dim(img)[3], 2)]

# ensure that image(s) have 3 dimenions
	if (length(dim(dapi)) == 2)	{
		dims <- dim(dapi)
		dim(dapi) <- c(dims, 1)
		dim(flr) <- c(dims, 1)
	}

# create nuclear mask, remove small and large nuclei
	xw <- nucMask(dapi, width = width, offset = offset, size = size,
			sigma = sigma, radius = radius, gamma = gamma)
	xw <- trimMask(xw, cutoff = cutoff, k = k)

# derive well descriptors
	well <-  well.info(basename(dirname(imageFiles)))$well[1]
	column <- well.info(well)$column
	row <- well.info(well)$row

# extract file name and number of frames
	fname <- paste(well, basename(imageFiles)[seq(2, length(imageFiles), 2)], sep = "/")
	nframes <- dim(xw)[3]

# assemble data
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
	return(list(data = res, mask = xw))
}

#
# process alternately stacked DAPI and fluorescent images in named file
# return data.frame as 'data' and nuclear mask as 'mask'
#
.processByStack <- function(imageStack, cutoff, k, width,
		offset, size, sigma, radius, gamma)
{
# read image files
	img <- suppressWarnings(readImage(imageStack))
	dapi <- img[,,seq(1, dim(img)[3], 2)]
	flr <- img[,,seq(2, dim(img)[3], 2)]

# ensure that image(s) have 3 dimenions
	if (length(dim(dapi)) == 2)	{
		dims <- dim(dapi)
		dim(dapi) <- c(dims, 1)
		dim(flr) <- c(dims, 1)
	}

# create nuclear mask, remove small and large nuclei
	xw <- nucMask(dapi, width = width, offset = offset, size = size,
			sigma = sigma, radius = radius, gamma = gamma)
	xw <- trimMask(xw, cutoff = cutoff, k = k)

# extract file name and number of frames
	fname <- basename(imageStack)
	nframes <- dim(xw)[3]

# assemble data
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
	return(list(data = res, mask = xw))
}
