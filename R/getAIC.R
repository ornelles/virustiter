#########################################################################################
#
# utility to evaluate AIC for fitted model(s) based on cut value(s)
#
#########################################################################################

getAIC <- function(df, cut, by = c("sequential", "column", "row")) {
	by <- match.arg(by)
	dfList <- lapply(cut, function(x) score(df, x))
	fmList <- lapply(dfList, function(x) getFit(x, by = by))
	if (class(class(fmList)[1])[1] == "list")
		ans <- lapply(fmList, function(x) sapply(x, AIC))
	else
		ans <- sapply(fmList, AIC)
	return(ans)
}
