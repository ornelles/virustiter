#########################################################################################
#
# utility to evaluate AIC for fitted model(s) based on cut value(s)
#
#########################################################################################
#' Obtain AIC for 'Cut' Value
#'
#' Apply a give \code{cut} value to create a fitted model (or models) and
#' return the Akaike's Criterion for the GLM fitted model. \code{getFit()}
#' is called for each value of \code{cut}.
#'
#' @param df Annotated \code{data.frame} or list of \code{data.frame}s with
#'   imaging results.
#' @param cut Numeric vector of cutoff values to test.
#' @param by A character string indicating the organization of the results as
#'   being a single result or by column or by rows.
#'
#' @return
#'
#' The AIC value associated with the fit for each value of cut or list of AIC
#' values if a list of data frames was provided.
#'
#' @export
#'
getAIC <- function(df, cut, by = c("sequential", "column", "row"))
{
	if (!is.data.frame(df))
		stop("'", deparse(substitute(df)), "' must be a data frame processed by score()")	
	by <- match.arg(by)
	dfList <- lapply(cut, function(x) score(df, x))
	fmList <- lapply(dfList, function(x) getFit(x, by = by))
	if (class(class(fmList)[1])[1] == "list")
		ans <- lapply(fmList, function(x) sapply(x, AIC))
	else
		ans <- sapply(fmList, AIC)
	return(ans)
}
