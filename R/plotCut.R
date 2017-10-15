#########################################################################################
# plotCut
# 
# Use lattice densityplot() to determine and display background values with find.bgnd()
#
# If cut is missing, it will be determined according to the 'by' parameter
#
# Note that this expects df to have been produced by readIJResults() If log==TRUE, log-
# transformed values are used identify breakpoints. Results are plotted and the
# selected cutoff values from get.bgnd() are returned as an invisible named vector if
# return.plot is FALSE, otherwise the lattice plot object is returned.
# 
#########################################################################################

plotCut <- function(df, cut, by = c("well", "row", "column", "file"), param = "val",
		smooth = 1, mult = 5, log = TRUE, main = NULL, as.table = TRUE,
		return.plot = FALSE, ...)
{
	if (missing(df)) {
		usage <- c("plotCut examples:",
			'  plotCut(df, by = "well", val, smooth=1, mult=5, log=TRUE)',
			'  plotCut(df) ## default values are same as above',
			'  plotCut(df, cut = 0.002)  ## uses cutoff value of 0.002',
			'  plotCut(df, groups=column, auto.key=T)')
		cat(usage, sep="\n")
		return(invisible(NULL))
	}

# parse arguments and perform error checking
	by <- match.arg(by)
	param <- gsub("\\\"","", deparse(substitute(param)))
	if (!param %in% names(df))
		stop("'", param, "' not in data set")
	if (!by %in% names(df))
		stop("'", by, "' not in data set")

	d.adj <- smooth	# to hand to density plot

	library(lattice)
	library(latticeExtra)

# create title
	if (is.null(main)) {
		main.text <- paste(deparse(substitute(df)), "  smooth=", signif(smooth,2),
			" mult=", signif(mult,2), sep="")
		main <- list(main.text, cex=1, font=1)
	}

# calculate background values and create strip labels 
	if (missing(cut))
		cut <- do.call(getCut, list(df, by, param, mult, log))
	else {
		labs <- as.character(unique(df[[by]]))
		cut <- rep(cut, length.out = length(labs))
		names(cut) <- labs
	}
	strip.labels <- paste(names(cut), signif(cut, 2), sep=" at ")
	form <- as.formula(paste("~", param, "|", as.factor(by)))
	xlist <- list()	# for log argument in scales
	if (log == TRUE) {
		xlist <- list(log = 10)
		cut <- log10(cut)
	}
	obj <- densityplot(form, df,
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
