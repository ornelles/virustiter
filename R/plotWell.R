#' Plot Results by Well with Lattice Graphics
#'
#' Plot a representation of the selected wells showing positive and
#' negative cells with \code{\link[lattice]{xyplot}}.
#'
#' @param df Annotated \code{data.frame} with imaging results.
#' @param myWell Character vector of wells to be plotted.
#' @param byFrame A \code{logical} value to separate frames.
#' @param cex Numeric value as factor by which to expand the plot symbol.
#' @param invert.y A \code{logical} value to invert y coordinates.
#' @param ... Additional arguments handed to \code{xyplot}.
#'
#' @return
#'
#' The plot is returned as an invisible \code{lattice} object.
#'
#' @import lattice
#'
#' @export
#'
plotWell <- function(df, myWell = levels(df$well), byFrame = TRUE,
	cex = 1, invert.y = TRUE, ...)
{
	if (!"well" %in% names(df))
		stop("requires variable 'well' in ", deparse(substitute(df)))
	if (!"positive" %in% names(df))
		df$positive <- FALSE
	if (invert.y)
		df$ym <- -df$ym

	myWell <- well.info(myWell)$well
	adj.cex <- cex
	myCex <- scale(sqrt(df$area), center = FALSE)
	leg.text <- as.character(round(quantile(df$area)))
	myCol <- trellis.par.get("superpose.symbol")$col
	myPch <- trellis.par.get("superpose.symbol")$pch
	if (byFrame == TRUE)
		form <- as.formula("ym ~ xm | well:factor(frame)")
	else
		form <- as.formula("ym ~ xm | well")
	obj <- xyplot(form, df, subset = well%in%myWell,
		groups = positive,
		cex = myCex*adj.cex,
		panel = function(x, y, ..., cex, subscripts, groups) {
					panel.xyplot(x, y, ..., cex = cex[subscripts],
					col = ifelse(groups[subscripts] == "FALSE", myCol[1], myCol[2]))
		},
		key = list(space = "right", adj = 1, between = 0.25,
					title = "Area (px)", cex.title = 1,
					text = list(leg.text),
					points = list(pch = myPch[1], cex = adj.cex*quantile(myCex), col = myCol[1]),
					points = list(pch = myPch[2], cex = adj.cex*quantile(myCex), col = myCol[2])),
		as.table = TRUE, aspect = "iso", xlab = "", ylab = "", scales = list(draw = FALSE),
		...)
	print(obj)
	invisible(obj)
}
