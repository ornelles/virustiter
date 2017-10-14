#########################################################################################
# getCut
#
# derived from generic code to determine cut (background) values for 'param' (val) with
# cut points from find.bgnd() by control wells (x==0 or type=="control"), or for each
# well, row, or column in the data.frame
#
# Dependencies
#   find.bgnd()	calculates likely cutoff point for background from density profile
#
#########################################################################################

getCut <- function(df, by=c("control","well","row","column"), param=val, mult=5, log=TRUE)
{
	if (missing(df)) {
		usage <- c("getCut examples:",
			'  getCut(df, "control", val, mult=5, log=TRUE)',
			'  getCut(df) # same as above',
			'  getCut(df, "row", mult=3, log=FALSE)')
		cat(usage, sep="\n")
		return(invisible(NULL))
	}
	if (!is.data.frame(df))
		stop(deparse(substitute(df)), " must be a data.frame")

# parse arguments and perform error checking
	by <- match.arg(by)
	param <- gsub("\\\"","", deparse(substitute(param)))
	if (!param %in% names(df))
		stop("'", param, "' not in data set")

# create local copy of relevant data
	if (by == "control") {				# special case of by "control"
		if (any(df$type=="control"))
			temp <- data.frame(g=TRUE, y=subset(df, type=="control")[[param]])
		else if (any(df$x==0))
			temp <- data.frame(g=TRUE, y=subset(df, x==0)[[param]])
		else
			stop('require type=="control" or x values of 0 to use by="control" option')
	}
	else {
		temp <- df[c(by, param)]
		names(temp) <- c("g", "y")
	}
	res <- aggregate(y ~ g, temp, function(v) find.bgnd(v, mult=mult, log=log))
	ret <- res$y
	names(ret) <- if(by=="control") "control" else res$g
	return(ret)
}
