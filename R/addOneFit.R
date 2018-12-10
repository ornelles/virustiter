#' Add Single Fitted Model to Plot
#'
#' Add a single fitted model plot to an existing plot made with base graphics
#'
#' Prepare a single diagnostic plot with more labeling options than \code{plotFit()}.
#'
#' @param fm Fitted model or list of fitted models from \code{getFit()}.
#' @param pch.col Pch color passed to \code{plot}.
#' @param line.col Color for best-fit line (appropriate for \code{par("col")}.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @details
#'
#' Base graphics are used to add a plot to an existing plot produced by 
#' \code{plotFit} or \code{plotOneFit}. The first argument is
#' a single fitted model from \code{getFit()}. The function calls
#' \code{getTiter()} to obtain the fit and confidence intervals.
#'
#' @return
#'
#' No value is returned. This function is called for the side-effect of 
#' adding to an existing plot with \code{points()} and \code{lines()}.
#'
#' @export
#'  
addOneFit <- function(fm, line.col=4, ref.col=4, pch.col=line.col, ...) {
	moi <- exp(fm$model[[2]])		# model data.frame holds values used for fit
	y <- prop.table(fm$model[[1]],1)[,1]
	res <- fm$data  # entire data.frame handed to glm()
	cf <- getTiter(fm)

	xlo <- min(moi[moi > 0])
	xhi <- max(moi)
	xp <- exp(seq(log(xlo), log(xhi), length=101))
	yp <- predict(fm, data.frame(x = xp), type="response")
	xpp <- cf[1]
	ypp <- 1 - exp(-1)

	points(y ~ moi, subset = moi > 0, col = pch.col, ...)
	lines(xp, yp, col = line.col)
	lines(c(xlo, xpp, xpp),c(ypp, ypp, -0.02), lty = 2, col = ref.col)
}
