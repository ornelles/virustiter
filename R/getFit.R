#########################################################################################
# getFit
#
# fit data in data.frame produced by tally() or by score() as sequential order, column, or row
# as indicated in argument 'by' and return glm model fit(s). Each fitted model has the
# unit assigned as an attribute with the text value
#
# by a two-column response in the binomial glm model, the fit is weighted by the total
# numbers of cases which is returned as prior.weights in the fitted model
#
#########################################################################################

getFit <- function(obj, by = c("sequential", "column", "row"))
{
	if (all(c("pos", "neg", "moi") %in% names(obj)))		# data.frame produced by tally()
		res <-  obj
	else if (all(c("positive", "well", "moi") %in% names(obj)))	# data.frame from score()
		res <- tally(obj)
	else
		stop(deparse(substitute(obj)), " does not to have been produced by tally() or score()")

	if ("unit" %in% names(res)) unit <- as.character(res$unit[1])
	else unit <- "none"
	by <- match.arg(by)

	fmList <- list()
	if (by == "sequential") {
		fm <- try(glm(cbind(pos, neg) ~ offset(log(moi)), data = res, subset = moi > 0,
				family=binomial("cloglog")))
		if (class(fm)[1] == "try-error")
			fmList[[1]] <- NULL
		else
			fmList[[1]] <- fm
	}
	else if (by == "column") {
		stopifnot("column" %in% names(res))
		for (i in levels(as.factor(res[[by]]))) {
			fm <- try(glm(cbind(pos, neg) ~ offset(log(moi)),
					family=binomial("cloglog"), data = res,
					subset = res[[by]] == i & moi > 0), silent=TRUE)
			if (class(fm)[1] == "try-error")
				fmList[[i]] <- NULL
			else
				fmList[[i]] <- fm
		}
	}
	else {	# by == "row"
		stopifnot("row" %in% names(res))
		for (i in levels(as.factor(res[[by]]))) {
			fm <- try(glm(cbind(pos, neg) ~ offset(log(moi)),
					family=binomial("cloglog"), data = res,
					subset = res[[by]] == i & moi > 0), silent=TRUE)
			if (class(fm)[1] == "try-error")
				fmList[[i]] <- NULL
			else
				fmList[[i]] <- fm
		}
	}
	for (i in seq_along(fmList))
		attr(fmList[[i]], "unit") <- unit 
	
# return fit or list of fits
	if (length(fmList)==1)
		return(fmList[[1]])
	else 
		return(fmList)
}

