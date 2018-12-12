#' Show Histogram of Mean Fluorescence Intensity with Background Cutoff
#'
#' Display a histogram of each well or file with \code{lattice} graphics showing
#' the selected background cutoff value or values.
#'
#' @param df Annotated \code{data.frame} with imaging results.
#' @param bgnd Numeric vector of background values. If missing, \code{getBgnd()}
#'    will be called with default parameters. 
#' @param by Character vector indicating grouping where "default" will use
#'   \code{well} if present or \code{file}. This value is also passed to
#'   \code{getBgnd} if necessary.
#' @param mult Numeric value passed to \code{getBgnd()} to scale bgnd.
#' @param log A \code{logical} value passed to \code{getBgnd()}.
#' @param param Name of the variable to be analyzed as a character string; 
#'   typically "mfi" or "val".
#' @param main Optional character string to serve as plot title.
#' @param as.table A \code{logical} value passed to \code{histogram()}.
#' @param layout An optional numeric vector to specify layout of histogram,
#'   passed to \code{histogram}.
#' @param ... Additional arguments handed to \code{histogram()}.
#'
#' @details
#'
#' This presents similar representation of the data as \code{plotDens())} and
#' can be used to examine the uniformity of results from an imaging experiment
#' and to iteratively check the \code{mult} argument provided to
#' \code{getBgnd()}. 
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
plotHist <- function(df, bgnd, by = c("default", "well", "file", "row", "column"),
		mult = 2.5, log = TRUE, param = "mfi", main = NULL, as.table = TRUE,
		layout = NULL, ...)
{
	if (missing(df)) {
		usage <- c("plotHist examples:",
			'  plotHist(df)      # calculates and plots default bgnd values in df',
			'  plotHist(df, bgnd) # where bgnd is explicitly provided')
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
			stop("'well' or 'file' is not in data set")
	}
	else if (!by %in% names(df))
		stop("'", by, "' is not in data set")

	if (!param %in% names(df))
		stop(deparse(substitute(param)), " is not in data set")

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

# create strip labels
	strip.labels <- paste(names(bgnd), signif(bgnd, 2), sep = " at ")

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
		bgnd <- log10(bgnd)
	}
	form <- as.formula(paste("~", param, "|", by))
	obj <- histogram(form, data = df, main = main,
		layout = layout, n = 64, as.table = as.table,
		panel = function(x, ..., subscripts)
		{
			panel.histogram(x,  ...)
			idx <- unique(df[[by]][subscripts])
			panel.abline(v = bgnd[idx], col = 2)
		},
		scales = list(x = xlist, y = list(relation = "free", rot = 0)),
		strip = strip.custom(factor.levels = strip.labels,
			par.strip.text = list(cex = 0.9)),
		xscale.components = xscale.components.log10ticks,
		...)
	plot(obj)
	invisible(obj)
}
