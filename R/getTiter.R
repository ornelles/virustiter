#########################################################################################
# getTiter
# extract titer with confidence intervals from fit (or list of fits)
#
#########################################################################################

getTiter <- function(fm)
{
	.coef <- function(fm) {
		cf <- numeric(3)
		if (is.null(fm))
			cf <- c(NA, NA, NA)
		else {
			cf[1] <- exp(-coef(fm))
			cfVal <- try(confint(fm))
			if(class(cfVal)[1] != "try-error")
				cf[c(3,2)] <- exp(-cfVal)
		}
		names(cf) <- c("est", "lo.CI", "hi.CI")
		cf
	}

	if ("list" %in% class(fm))
		ret <- t(sapply(fm, .coef))
	else
		ret <- .coef(fm)
	return(ret)
}

