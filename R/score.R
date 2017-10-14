#########################################################################################
# score
#
# Assign positive values to those exceeding the cutoff values in 'cut'
# 'cut' can be a single value or named vector for well, row, or column from getCut()
# if 'cut' is missing, default settings (control) of getCut() will be used 
#
# Arguments
#	df			data frame from readIJResults()
#	cut			cutoff value(s)
#
# Returns original data frame with 'positive' scored as TRUE or FALSE
#
#########################################################################################

score <- function(df, cut)
{
	if (missing(df)) {
		usage <- c("score examples:",
			'  score(df, cut)        ## cut holds cutoff values',
			'  score(df, getCut(df)) ## cutoff values obtained from control wells',
			'  score(df, 256)        ## single value used as cutoff')
		cat(usage, sep="\n")
		return(invisible(NULL))
	}
	stopifnot(c("val") %in% names(df))
	if (missing(cut)) {
		message("using default parameters in getCut() for 'cut' parameter")
		cut <- getCut(df)
	}

	if (length(cut) == 1)			# single values
		df$positive <- df$val > cut
	else if (all(names(cut) %in% levels(df$well)))
		df$positive <- df$val > cut[df$well]
	else if (all(names(cut) %in% levels(df$row)))
		df$positive <- df$val > cut[df$row]
	else if (all(names(cut) %in% levels(df$column)))
		df$positive <- df$val > cut[df$column]
	else
		stop("cut must be a single value or a named vector (well, row, or column)")
	return(df)
}
