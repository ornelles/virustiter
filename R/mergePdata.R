#########################################################################################
# mergePdata
#
# Add phenotype data to image data obtained by parseImages()
# phenoData must have "moi", either "well" or "file"
#
#########################################################################################

mergePdata <- function(phenoData, imageData)
{
# determine data type, check arguments and harmonize well names
	stopifnot("moi" %in% names(phenoData))

# for separate images in folders (multi-well)
	if ("well" %in% names(imageData)) {
		stopifnot("well" %in% names(phenoData))
		phenoData$well <- factor(well.info(phenoData$well)$well) # harmonize
		stopifnot(levels(imageData$well) %in% levels(phenoData$well))
	# add row and column information to phenotype data
		phenoData$column <- well.info(phenoData$well)$column
		phenoData$row <- well.info(phenoData$well)$row
	}
	else
		stopifnot("file" %in% names(phenoData))

# check value for unit and assign control wells
	if (!"unit" %in% names(phenoData))
		phenoData$unit <- "unspecified"
	type <- rep("standard", nrow(phenoData))
	type[phenoData$moi == 0] <- "control"
	phenoData$type <- type

# remove any variables in imageData that are present in phenoData EXCEPT
# for variables used to merge data frames: 'file' and/or 'well'
	vars <- names(imageData)[!names(imageData) %in% names(phenoData)]
	if ("file" %in% names(imageData)) vars <- c(vars, "file")
	if ("well" %in% names(imageData)) vars <- c(vars, "well")
	vars <- unique(vars)
	imageData <- imageData[vars]
	res <- merge(phenoData, imageData)

# convert strings to factors
	sel <- sapply(res, is.character)
	res[sel] <- lapply(res[sel], as.factor)

# reorganize data
	pdnames <- c("directory","file","column","row","well","frame","type","unit")
	first <- pdnames[pdnames %in% names(res)]
	last <- names(res)[!names(res) %in% pdnames]
	res <- res[c(first, last)]

# exclude unused levels (allows pd to hold information than required)
	res <- droplevels(res)
	rownames(res) <- NULL
	return(res)
}
