#' Add Single Fitted Model to Plot
#'
#' Add a single fitted model plot to an existing plot made with base graphics
#'
#' @param fm single fitted model from \code{\link{getFit}}.
#' @param pch,col,col.pch plot character and color for data points and the 
#'  default color for the fitted line. Set \code{pch} to \code{NA}
#'  to exclude plotting characters.
#' @param lty.fit,col.fit line type and color for the GLM best-fit line.  Set
#'  \code{lty.fit} to \code{NA} to exclude the best-fit line.
#' @param lty.ref,col.ref line type and color for the value on the x-axis
#'  that intersect the 63% value on the best-fit line. Set \code{lty.ref} to
#'  \code{NA} to exclude the reference line.
#' @param ... Additional arguments passed to \code{\link{plot}}.
#'
#' @details
#'
#' Base graphics are used to add a plot to an existing plot produced by 
#' \code{\link{plotFit}}. The first argument is
#' a single fitted model from \code{\link{getFit}}. The function calls
#' \code{\link{getTiter}} to obtain the fit and confidence intervals.
#'
#' @return
#'
#' No value is returned. This function is called for the side-effect of 
#' adding to an existing plot with \code{\link{points}} and \code{\link{lines}}.
#'
#' @export
#'  
addOneFit <- function(fm, pch = 1, col = 2, col.pch = col, lty.fit = 1,
	col.fit = col, lty.ref = 2, col.ref = "gray",  ...)
{
	# bookkeeping 
		x <- exp(fm$model[[2]])		# model data.frame holds values used for fit
		y <- prop.table(fm$model[[1]],1)[,1]
		res <- fm$data  # entire data.frame handed to glm()
	
		xlo <- min(x[x > 0])
		xhi <- max(x)
		xp <- exp(seq(log(xlo), log(xhi), length = 51))
		yp <- predict(fm, data.frame(x = xp), type = "response")
		xpp <- getTiter(fm, NULL)
		ypp <- 1 - exp(-1)
	
	# points
		if (!is.na(pch))
			points(y ~ x, subset = x > 0, pch = pch, col = col.pch, ...)
	# fitted line
		if (!is.na(lty.fit))
			lines(xp, yp, col = col.fit, lty = lty.fit, ...)
	# reference line
		if (!is.na(lty.ref))
			lines(c(xlo, xpp, xpp), c(ypp, ypp, -0.02), lty = lty.ref,
				col = col.ref, ...)
}
