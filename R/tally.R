#########################################################################################
# tally
#
# tally positive and negative values data.frame from readIJResults() by well
# merge with optional phenotype data frame, pd
#
# returns "result" data.frame with dname, well, moi, pos, neg and y and unit
#
#########################################################################################

tally <- function(df, pd)
{
	stopifnot(c("well", "positive", "moi") %in% names(df))

# extract data frame name

	if (!is.null(df$dname) & nlevels(df$dname)==1)
		dname <- levels(df$dname)
	else
		dname <- "unknown"

# tally positive and create results data.frame

	pos <- tapply(df$positive==TRUE, df$well, sum)
	neg <- tapply(df$positive==FALSE, df$well, sum)
	y <- pos/(pos + neg)
	moi <- sapply(names(pos),function(v) df$moi[df$well==v][1])
	unit <- df$unit[1]
	well <- names(pos)
	row <- well.info(well)$row
	column <- well.info(well)$column
	
	
	res <- data.frame(dname, well, row, column, moi, unit, pos, neg, y)
	if (!missing(pd)) {
		stopifnot("well" %in% names(pd))
		res <- merge(pd, res)
	}
	return(res)
}

