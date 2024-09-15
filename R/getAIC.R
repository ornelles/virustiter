
#' Obtain AIC for Given Background Value
#'
#' Apply a given \code{bgnd} value to create a fitted model (or models) and
#' return the Akaike's Criterion for the GLM fitted model. \code{\link{getFit}}
#' is called for each value of \code{bgnd}.
#'
#' @param df annotated \code{data.frame} or list of \code{data.frame}s with
#'   imaging results.
#' @param bgnd numeric vector of background values to test.
#' @param by (optional) name of the variable in the results \code{data.frame}
#'   that will be used to split the results before fitting
#'
#' @return
#'
#' The AIC value associated with the fit for each value of bgnd or list of AIC
#' values if a list of data frames was provided.
#'
#' @export
#'
getAIC <- function(df, bgnd, by)
{
	if (!is.data.frame(df))
		stop("'", deparse(substitute(df)),
			"' must be a data frame processed by score()")		
	if (missing(by))
		by <- NULL
	else if (by %in% names(df))
		by <- by
	else
		stop("unable to find ", deparse(substitute(by)), " in data.frame")

	dfList <- lapply(bgnd, function(x) score(df, x))
	fmList <- lapply(dfList, function(x) getFit(x, by = by))
	if (class(class(fmList)[1])[1] == "list")
		ans <- lapply(fmList, function(x) sapply(x, AIC))
	else
		ans <- sapply(fmList, AIC)
	return(ans)
}
