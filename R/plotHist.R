#' Show Histogram of Mean Fluorescence Intensity with Background Cutoff
#'
#' Display a histogram of each well or file with \code{lattice} graphics showing
#' the selected background cutoff values.
#'
#' @param df Annotated \code{data.frame} with imaging results.
#' @param bgnd Numeric vector of background values. If missing, \code{\link{getBgnd}}
#'    will be called with parameters \code{by, param, mult,} and \code{log}.
#' @param param Character string identifying the variable to be analyzed. Also
#'   passed to \code{\link{getBgnd}} if required.
#' @param panel Optional character string defining the \code{lattice} panels,
#'   typically \code{"well"} or \code{"file"}. 
#' @param log Optional \code{logical} or \code{numeric} value to transform
#'   \code{'param'} values. Also passed to \code{\link{getBgnd}} if required. 
#' @param by,mult Additional parameters passed to \code{\link{getBgnd}} if required.
#' @param main Optional character string to serve as plot title. If \code{NULL},
#'   the system date will be used.
#' @param as.table A \code{logical} value passed to \code{histogram()}.
#' @param layout An optional numeric vector to specify layout of histogram,
#'   passed to \code{\link[lattice]{histogram}}.
#' @param ... Additional arguments passed to \code{\link[lattice]{histogram}}.
#'
#' @details
#'
#' A histogram of image intensity in \code{'param'} is plotted
#' with the selected background cutoff similar to the function
#' \code{\link{plotDens}}. Both functions can be used to examine the uniformity
#' of results from an imaging experiment and to interactively check the
#' paramaters passed to \code{\link{getBgnd}} to determine a suitable background
#' cutoff value. 
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
plotHist <- function(df, bgnd, param = "mfi", panel, log = TRUE, by = NULL, 
		mult = NULL, main = NULL, as.table = TRUE, layout = NULL, ...)
{
	if (missing(df)) {
		usage <- c("plotHist examples:",
			'  plotHist(df)       # calculate and plot default bgnd values in df',
			'  plotHist(df, bgnd) # where bgnd is explicitly provided')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}
	requireNamespace("lattice", quietly = TRUE)
	requireNamespace("latticeExtra", quietly = TRUE)

# check 'param' argument
	if (!param %in% names(df))
		stop(deparse(substitute(param)), " is not in data set")

# process 'bgnd' argument
	if (missing(bgnd)) {
		argNames <- names(formals("getBgnd"))
		sel <- sapply(mget(argNames), is.null)
		bgnd <- do.call("getBgnd", mget(argNames)[!sel])
	}

# process 'panel' argument
	if (missing(panel)) {
		if ("well" %in% names(df)) panel <- "well"
		if ("well" %in% names(df)) panel <- "well"
		else if ("file" %in% names(df)) panel <- "file"
		else stop("unable to assign panel group as 'well' or 'file'")
	}
	else if (!is.character(panel))
		stop("'panel' must be a character vector")
	else if (!panel %in% names(df))
		stop("'", panel, "' not in data set")

# determine factor index for strip labels from 'bgnd'
	if (is.null(names(bgnd)))
		index <- NA
	else if(length(names(bgnd)) == 1 && names(bgnd) == "control")
		index <- NA
	else { # search for match among factors
		sel <- sapply(lapply(df, levels), function(v) all(names(bgnd) %in% v))
		index <- names(which(sel)[1]) # first one that matches
	}

# assign names to 'bgnd' to use as strip labels
	if (!is.na(index)) {
		mat <- unique(df[c(panel, index)]) # two column matrix
		ord <- order(mat[[index]]) # preserve order of factors
		mat <- mat[ord,]
		bgnd <- bgnd[as.character(mat[[2]])]
		lab.panel <- as.character(mat[[1]])
		lab.bgnd <- names(bgnd)
		if (all(lab.bgnd %in% lab.panel))
			names(bgnd) <- lab.bgnd
		else
			names(bgnd) <- paste(lab.panel, lab.bgnd)
	}
	else { # single background value provided
		ord <- order(levels(as.factor(df[[panel]]))) # preserve order of factors
		lab.panel <- as.character(unique(df[[panel]]))[ord]
		bgnd <- rep(bgnd, length(lab.panel))
		names(bgnd) <- lab.panel
	}

# create strip labels
	strip.labels <- paste(names(bgnd), signif(bgnd, 2), sep = " at ")

# assemble lattice plot
	if (is.null(layout))
		layout <- c(1, nlevels(df[[panel]]))
	if (is.null(main)) {
		main.text <- Sys.Date()
		main <- list(main.text, cex = 1, font = 1)
	}

# process 'log' argument
	if (identical(log, NULL)) logsc <- FALSE
	else if (identical(log, FALSE)) logsc <- FALSE
	else if (identical(log, TRUE)) logsc <- 10
	else if (log == 1) logsc <- 10
	else logsc <- as.numeric(log)

# prepare x scale
	xlist <- list()	
	if (logsc != 0) {
		xlist <- list(log = logsc)
		bgnd <- log(bgnd, logsc)
	}
	if (logsc == 10)
		xsc <- latticeExtra::xscale.components.log10ticks
	else if (is.numeric(logsc))
		xsc <- latticeExtra::xscale.components.logpower
	else
		xsc <- latticeExtra::xscale.components.default

# create lattice formula and object
	form <- as.formula(paste("~", param, "|", panel))

	obj <- histogram(form, data = df, main = main, layout = layout, n = 64,
		panel = function(x, ..., subscripts) {
			panel.histogram(x,  ...)
			index <- unique(df[[panel]][subscripts])
			panel.abline(v = bgnd[index], col = 2)
		},
		scales = list(x = xlist, y = list(relation = "free", rot = 0)),
		xscale.components = xsc, as.table = as.table,
		strip = strip.custom(factor.levels = strip.labels,
			par.strip.text = list(cex = 0.9)),
		...)

	plot(obj)
	invisible(obj)
}
