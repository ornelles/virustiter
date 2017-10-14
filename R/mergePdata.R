#########################################################################################
# mergePdata
#
# Add phenotype data to image data obtained by parseImages()
# phenoData must have at least "well", "x" (moi), and should have "unit"
#
#########################################################################################

mergePdata <- function(phenoData, imageData) {
	stopifnot("well" %in% names(phenoData))
	stopifnot("x" %in% names(phenoData))
	stopifnot("well" %in% names(imageData))

# harmonize phenoData well names, add row and column information, clean up
	phenoData$well <- factor(well.info(phenoData$well)$well)
	stopifnot(all(levels(imageData$well) %in% levels(phenoData$well)))
	phenoData$column <- well.info(phenoData$well)$column
	phenoData$row <- well.info(phenoData$well)$row

	if (!"unit" %in% names(phenoData))
		phenoData$unit <- "unknown"
	type <- rep("standard", nrow(phenoData))
	type[phenoData$x==0] <- "control"
	phenoData$type <- type
	
	sel <- names(phenoData) %in% c("well","column","row", "unit", "type")
	phenoData <- cbind(phenoData[sel],phenoData[!sel])

	res <- merge(phenoData, imageData)
	sel <- sapply(res, is.character)
	res[sel] <- lapply(res[sel], as.factor)
	rownames(res) <- NULL
	return(res)
}
