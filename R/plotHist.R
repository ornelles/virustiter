#########################################################################################
# plotHist
#
# display histogram for each well or file with lattice graphics with background cutoff
# 'cut' value as a single value or named vector from getCut()
#
# if missing and 'positive' is present, the maximum value will be shown
#
#########################################################################################
#' Show Cutoff with Histogram of Mean Fluorescence Intensity
#'
#' Display a histogram of each well or file with \code{lattice} graphics showing
#' the selected background cutoff value or values.
#'
#' @param df Annotated \code{data.frame} with imaging results.
#' @param cut Numeric vector of cutoff values. If missing, \code{getCut} will be
#'   called with default parameters. 
#' @param by Character vector indicating grouping where "default" will use
#'   \code{well} if present or \code{file}. This value is also passed to
#'   \code{getCut} if necessary.
#' @param mult Numeric value passed to \code{getCut} to scale cutoff.
#' @param log A \code{logical} value passed to \code{getCut}.
#' @param param Name of the variable to be analyzed as a character string; 
#'   typically "mfi" or "val".
#' @param main Optional character string to serve as plot title.
#' @param as.table A \code{logical} value passed to \code{histogram}.
#' @param layout An optional numeric vector to specify layout of histogram,
#'   passed to \code{histogram}.
#' @param ... Additional arguments handed to \code{histogram}.
#'
#' @details
#'
#' This presents similar representation of the data as \code{plotDensity} and
#' can be used to examine the uniformity of results from an imaging experiment
#' and to iteratively check the \code{mult} argument provided to \code{getCut}. 
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
plotHist <- function(df, cut, by = c("default", "well", "file", "row", "column"),
		mult = 5, log = TRUE, param = "mfi", main = NULL, as.table = TRUE,
		layout = NULL, ...)
{
	if (missing(df)) {
		usage <- c("plotHist examples:",
			'  plotHist(df)      # calculates and plots default cut values in df',
			'  plotHist(df, cut) # where cut is explicitly provided')
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

# calculate background cutoff value and assign names 
	if (missing(cut))
		cut <- do.call(getCut, list(df, by, param, mult, log))
	else {
		labs <- as.character(unique(df[[by]]))
		cut <- rep(cut, length.out = length(labs))
		names(cut) <- labs
	}

# create strip labels
	strip.labels <- paste(names(cut), signif(cut, 2), sep = " at ")

# assemble lattice plot
	if (is.null(layout))
		layout <- c(1, nlevels(df[[by]]))
	if (is.null(main)) {
		main.text <- paste0(deparse(substitute(df)), " mult = ", signif(mult, 2))
		main <- list(main.text, cex = 1, font = 1)
	}
	xlist <- list()	# for log argument in scales
	if (log == TRUE) {
		xlist <- list(log = 10)
		cut <- log10(cut)
	}
	form <- as.formula(paste("~", param, "|", by))
	obj <- histogram(form, data = df, main = main,
		layout = layout, n = 64, as.table = as.table,
		panel = function(x, ..., subscripts)
		{
			panel.histogram(x,  ...)
			idx <- unique(df[[by]][subscripts])
			panel.abline(v = cut[idx], col = 2)
		},
		scales = list(x = xlist, y = list(relation = "free", rot = 0)),
		strip = strip.custom(factor.levels = strip.labels,
			par.strip.text = list(cex = 0.9)),
		xscale.components = xscale.components.log10ticks,
		...)
	plot(obj)
	invisible(obj)
}
