#' Show Density Plot of Mean Fluorescence Intensity with Background Cutoff
#'
#' Display a \code{densityplot} of each well or file with \code{lattice} graphics
#' showing the selected background cutoff value or values.
#'
#' @param df Annotated \code{data.frame} with imaging results.
#' @param bgnd Numeric vector of bgnd values. If missing, \code{getBgnd()}
#'   will be called with default parameters. 
#' @param by Character vector indicating grouping where "default" will use
#'   \code{well} if present or \code{file}. This value is also passed to
#'   \code{getBgnd()} if necessary.
#' @param smooth Numeric value passed to the \code{density()} function.
#' @param mult Numeric value passed to \code{getBgnd()} to scale bgnd.
#' @param log A \code{logical} value passed to \code{getBgnd()}.
#' @param main Optional character string to serve as plot title.
#' @param as.table A \code{logical} value passed to \code{histogram()}.
#' @param param Name of the variable to be analyzed as a character string; 
#'   typically "mfi" or "val".
#' @param ... Additional arguments handed to \code{densityplot()}.
#'
#' @details
#'
#' This presents representation of the data like \code{plotHist()} and can
#' be used to examine the uniformity of results from an imaging experiment and
#' to iteratively check the \code{mult} argument provided to \code{getBgnd()}.
#' The background values are incorporated into the strip labels. 
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
plotDens <- function(df, bgnd, by = c("default", "well", "file", "row", "column"),
		smooth = 1, mult = 2.5, log = TRUE, main = NULL, as.table = TRUE,
		param = "mfi", ...)
{
	if (missing(df)) {
		usage <- c("plotDens examples:",
			'  plotDens(df, by = "well", smooth = 1, mult = 2, log = TRUE)',
			'  plotDens(df) ## default values are same as above',
			'  plotDens(df, bgnd = 0.002)  ## uses bgnd value of 0.002',
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

# calculate background cutoff value if necessary
	if (missing(bgnd))
		bgnd <- do.call("getBgnd", list(df, by, param, mult, log))

# assign names to bgnd to use as strip labels
	if (is.null(names(bgnd)))
		idx <- NA
	else {
		sel <- sapply(lapply(df, levels), function(v) all(names(bgnd) %in% v))
		idx <- names(which(sel)[1]) # first one that matches
	}
	if (!is.na(idx)) {
		mat <- unique(df[c(by, idx)]) # two column matrix
		bgnd <- bgnd[as.character(mat[[2]])]
		lab1 <- as.character(mat[[1]])
		lab2 <- names(bgnd)
		if (all(lab2 %in% lab1))
			names(bgnd) <- lab2
		else
			names(bgnd) <- paste(lab1, lab2)
	}
	else { # single background value provided
		labs <- as.character(unique(df[[by]]))
		bgnd <- rep(bgnd, length(labs))
		names(bgnd) <- labs
	}

# create plot title
	if (is.null(main)) {
		main.text <- paste(deparse(substitute(df)), "  smooth = ", signif(smooth, 2),
			" mult = ", signif(mult, 2), sep = "")
		main <- list(main.text, cex = 1, font = 1)
	}

# adjust strip labels to show cutoff values
	strip.labels <- paste(names(bgnd), signif(bgnd, 2), sep = " at ")
	form <- as.formula(paste("~", param, "|", by))
	xlist <- list()	# for log argument in scales
	if (log == TRUE) {
		xlist <- list(log = 10)
		bgnd <- log10(bgnd)
	}
	obj <- densityplot(form, data = df,
		scales = list(x = xlist, y = list(draw = FALSE, relation = "free")),
		main = main,
		panel = function(x, myBgnd = bgnd) {
			panel.densityplot(x, plot.points = FALSE)
			panel.abline(v = myBgnd[panel.number()], col = 2)},
		as.table = as.table,
		strip = strip.custom(factor.levels = strip.labels,
			par.strip.text = list(cex = 0.9)),
		xscale.components = xscale.components.log10ticks, ...)
	plot(obj)
	if (log == TRUE)	# return to linear scale
		bgnd <- 10^bgnd
	return(invisible(obj))	# return lattice plot
}
