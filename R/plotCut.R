#########################################################################################
# plotCut
# 
# Use lattice densityplot() to determine and display background values with findBgnd()
#
# If cut is missing, it will be determined according to the 'by' parameter
#
# Note that this expects df to have been produced by readIJResults() If log==TRUE, log-
# transformed values are used identify breakpoints. Results are plotted and the
# selected cutoff values from get.bgnd() are returned as an invisible named vector if
# return.plot is FALSE, otherwise the lattice plot object is returned.
# 
#########################################################################################

plotCut <- function(df, cut, by = c("default", "well", "file", "row", "column"),
		smooth = 1, mult = 5, log = TRUE, main = NULL, as.table = TRUE,
		return.plot = FALSE, ...)
{
	if (missing(df)) {
		usage <- c("plotCut examples:",
			'  plotCut(df, by = "well", smooth=1, mult=5, log=TRUE)',
			'  plotCut(df) ## default values are same as above',
			'  plotCut(df, cut = 0.002)  ## uses cutoff value of 0.002',
			'  plotCut(df, groups = column, auto.key = T)')
		cat(usage, sep="\n")
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
	d.adj <- smooth	# to hand to density plot

# calculate background cutoff value and create strip labels 
	if (missing(cut))
		cut <- do.call(getCut, list(df, by, "val", mult, log))
	else {
		labs <- as.character(unique(df[[by]]))
		cut <- rep(cut, length.out = length(labs))
		names(cut) <- labs
	}

# create plot title
	if (is.null(main)) {
		main.text <- paste(deparse(substitute(df)), "  smooth=", signif(smooth,2),
			" mult=", signif(mult,2), sep="")
		main <- list(main.text, cex=1, font=1)
	}
	strip.labels <- paste(names(cut), signif(cut, 2), sep=" at ")
	form <- as.formula(paste("~ val |", as.factor(by)))
	xlist <- list()	# for log argument in scales
	if (log == TRUE) {
		xlist <- list(log = 10)
		cut <- log10(cut)
	}
	obj <- densityplot(form, data = df,
		scales = list(x=xlist, y=list(draw=FALSE, relation="free")), main = main,
		panel = function(x, bgnd) panel.bgnd(x, bgnd = cut), as.table=as.table,
		strip=strip.custom(factor.levels=strip.labels, par.strip.text=list(cex=0.9)),
		xscale.components=xscale.components.log10ticks, ...)
	plot(obj)
	if (log == TRUE)	# return to linear scale
		cut <- 10^cut
	if (return.plot==FALSE)
		return(invisible(cut))	# return background values
	else
		return(invisible(obj))	# return lattice plot
}

#########################################################################################
#
# lattice panel function to plot density and cutoff value
#
#########################################################################################
##
## REVISE TO RESPECT CUT from PLOT CUT
##
panel.bgnd <- function(x, bgnd)
{
	panel.densityplot(x, plot.points = FALSE)
	panel.abline(v = bgnd[panel.number()], col = 2)
}
