#########################################################################################
# mergePdata
#
# Add phenotype data to image data obtained by parseImages()
# phenoData must have at least "well", "moi", and should have "unit"
#
#########################################################################################

mergePdata <- function(phenoData, imageData)
{
# check arguments and harmonize well names
	stopifnot("well" %in% names(imageData))
	stopifnot("well" %in% names(phenoData))
	stopifnot("moi" %in% names(phenoData))
	phenoData$well <- factor(well.info(phenoData$well)$well) # harmonize
	stopifnot(levels(imageData$well) %in% levels(phenoData$well))

# add row and column information to phenotype data
	phenoData$column <- well.info(phenoData$well)$column
	phenoData$row <- well.info(phenoData$well)$row

# check value for unit and assign control wells
	if (!"unit" %in% names(phenoData))
		phenoData$unit <- "unspecified"
	type <- rep("standard", nrow(phenoData))
	type[phenoData$moi == 0] <- "control"
	phenoData$type <- type
	
# reorganize phenotype data
	sel <- names(phenoData) %in% c("well", "column", "row", "unit", "type")
	phenoData <- cbind(phenoData[sel], phenoData[!sel])

# remove variables in imageData that are already present in phenoData
	keep <- names(imageData)[!names(imageData) %in% names(phenoData)]
	imageData <- imageData[c("well", keep)]
	res <- merge(phenoData, imageData)

# convert strings to factors
	sel <- sapply(res, is.character)
	res[sel] <- lapply(res[sel], as.factor)

# exclude unused levels (allows pd to hold information than required)
	res <- droplevels(res)
	rownames(res) <- NULL
	return(res)
}
