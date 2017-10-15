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
	
# reorganize phenotype data
	pdnames <- c("directory","file","well","frame","column","row","unit","type")
	sel <- names(phenoData) %in% pdnames
	phenoData <- cbind(phenoData[sel], phenoData[!sel])

# remove variables in imageData that are already present in phenoData
	keep <- names(imageData)[!names(imageData) %in% names(phenoData)]
	special <- c("file", "well")
	special <- special[special %in% names(imageData)]
	imageData <- imageData[c(special, keep)]
	res <- merge(phenoData, imageData)

# convert strings to factors
	sel <- sapply(res, is.character)
	res[sel] <- lapply(res[sel], as.factor)

# exclude unused levels (allows pd to hold information than required)
	res <- droplevels(res)
	rownames(res) <- NULL
	return(res)
}
