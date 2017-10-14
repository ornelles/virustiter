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
	stopifnot("well" %in% names(phenoData))
	phenoData$well <- factor(well.info(phenoData$well)$well)
	stopifnot(levels(imageData$well) %in% levels(phenoData$well))
	stopifnot("moi" %in% names(phenoData))
	stopifnot("well" %in% names(imageData))

# add row and column information to data
	phenoData$column <- well.info(phenoData$well)$column
	phenoData$row <- well.info(phenoData$well)$row

	if (!"unit" %in% names(phenoData))
		phenoData$unit <- "unknown"
	type <- rep("standard", nrow(phenoData))
	type[phenoData$moi == 0] <- "control"
	phenoData$type <- type
	
	sel <- names(phenoData) %in% c("well","column","row", "unit", "type")
	phenoData <- cbind(phenoData[sel], phenoData[!sel])

	res <- merge(phenoData, imageData)
	sel <- sapply(res, is.character)
	res[sel] <- lapply(res[sel], as.factor)
# exclude unused levels to allow pd to hold more than required
	res <- droplevels(res)
	rownames(res) <- NULL
	return(res)
}
