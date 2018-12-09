#' Show Cutoff Values on a Densityplot of Mean Fluorescence Intensity
#'
#' Display a \code{densityplot} of each well or file with \code{lattice} graphics
#' showing the selected background cutoff value or values.
#'
#' @param df Annotated \code{data.frame} with imaging results.
#' @param cut Numeric vector of cutoff values. If missing, \code{getCut} will be
#'   called with default parameters. 
#' @param by Character vector indicating grouping where "default" will use
#'   \code{well} if present or \code{file}. This value is also passed to
#'   \code{getCut} if necessary.
#' @param smooth Numeric value passed to the \code{density} function.
#' @param mult Numeric value passed to \code{getCut} to scale cutoff.
#' @param log A \code{logical} value passed to \code{getCut}.
#' @param main Optional character string to serve as plot title.
#' @param as.table A \code{logical} value passed to \code{histogram}.
#' @param param Name of the variable to be analyzed as a character string; 
#'   typically "mfi" or "val".
#' @param ... Additional arguments handed to \code{histogram}.
#'
#' @details
#'
#' This presents similar representation of the data as \code{plotHist} and can be
#' used to examine the uniformity of results from an imaging experiment and to
#' iteratively check the \code{mult} argument provided to \code{getCut}. The
#' cutoff values are incorporated into the strip labels. 
#'
#' @return
#'
#' The plot is returned as an invisible \code{lattice} object.
#'
#' @import
#' lattice
#' latticeExtra
#'
#' @export
#'  
plotDens <- function(df, cut, by = c("default", "well", "file", "row", "column"),
		smooth = 1, mult = 2.5, log = TRUE, main = NULL, as.table = TRUE,
		param = "mfi", return.plot = FALSE, ...)
{
	if (missing(df)) {
		usage <- c("plotDens examples:",
			'  plotDens(df, by = "well", smooth = 1, mult = 2, log = TRUE)',
			'  plotDens(df) ## default values are same as above',
			'  plotDens(df, cut = 0.002)  ## uses cutoff value of 0.002',
			'  plotDens(df, groups = column, auto.key = T)')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}
	library(lattice)
	library(latticeExtra)

# parse arguments and perform error checking
	by <- match.arg(by)
	if (by == "default") {
		if ("well" %in% names(df))
			by <- "well"
		else if ("file" %in% names(df))
			by <- "file"
		else
			stop("'well' and 'file' not in data set")
	}
	else if (!by %in% names(df))
		stop("'", by, "' not in data set")
	if (!param %in% names(df))
		stop(deparse(substitute(param)), " not in data set")
	d.adj <- smooth	# to hand to density plot

# calculate background cutoff value and create strip labels 
	if (missing(cut))
		cut <- do.call(getCut, list(df, by, param, mult, log))
	else {
		labs <- as.character(levels(df[[by]]))
		cut <- rep(cut, length.out = length(labs))
		names(cut) <- labs
	}

# create plot title
	if (is.null(main)) {
		main.text <- paste(deparse(substitute(df)), "  smooth = ", signif(smooth,2),
			" mult = ", signif(mult, 2), sep = "")
		main <- list(main.text, cex = 1, font = 1)
	}

# adjust strip labels to show cutoff
	strip.labels <- paste(names(cut), signif(cut, 2), sep = " at ")
	form <- as.formula(paste("~", param, "|", by))
	xlist <- list()	# for log argument in scales
	if (log == TRUE) {
		xlist <- list(log = 10)
		cut <- log10(cut)
	}
	obj <- densityplot(form, data = df,
		scales = list(x = xlist, y = list(draw = FALSE, relation = "free")),
		main = main,
		panel = function(x, bgnd = cut) {
			panel.densityplot(x, plot.points = FALSE)
			panel.abline(v = bgnd[panel.number()], col = 2)},
		as.table = as.table,
		strip = strip.custom(factor.levels = strip.labels,
			par.strip.text = list(cex = 0.9)),
		xscale.components = xscale.components.log10ticks, ...)
	plot(obj)
	if (log == TRUE)	# return to linear scale
		cut <- 10^cut
	return(invisible(obj))	# return lattice plot
}
