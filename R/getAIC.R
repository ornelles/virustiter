#########################################################################################
#
# utility to evaluate AIC for fitted models based on cut value
#
#########################################################################################

getAIC <- function(df, cut, by = c("sequential","column","row")) {
	by <- match.arg(by)
	df <- score(df, cut)
	fm <- getFit(df, by = by)
	if (class(fm)[1] == "list")
		return(sapply(fm, AIC))
	else
		return(AIC(fm))
}
