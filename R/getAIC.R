#########################################################################################
#
# utility to evaluate AIC for fitted models based on cut value
#
#########################################################################################

getAIC <- function(df, cut, by="sequential") {
	df <- score(df, cut)
	fm <- getFit(df, by=by)
	if (class(fm)[1]=="list")
		return(sapply(fm, AIC))
	else
		return(AIC(fm))
}
