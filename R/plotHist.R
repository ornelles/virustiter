#########################################################################################
# plotHist
#
# display histogram for each well or file with lattice graphics with background cutoff
# 'cut' value as a single value or named vector from getCut()
#
#########################################################################################

plotHist <- function(df, cut = NULL, layout = NULL, ...) {
	if (missing(df)) {
		usage <- c("plotHist examples:",
			'  plotHist(df)      ## uses "positive" values in df',
			'  plotHist(df, cut) ## where cut was from getCut()')
		cat(usage, sep="\n")
		return(invisible(NULL))
	}
	library(lattice)
# select well or file as grouping variable
	if ("well" %in% names(df))
		group <- df$well
	else
		group <- df$file

	if (is.null(cut))					# use existing 'positive' assignment
		cut.points <- with(df, tapply(val, list(positive, group), max))[1,]
	else if (length(cut) == 1)			# single value
		cut.points <- rep(cut, nlevels(group))
	else if (all(names(cut) %in% levels(group)))
		cut.points <- cut
	else if (all(names(cut) %in% levels(df$row)))
		cut.points <- rep(cut, each=nlevels(df$column))
	else if (all(names(cut) %in% levels(df$column)))
		cut.points <- rep(cut, nlevels(df$row))
	else
		stop("cut must be null, a single number, or a named vector (file, well, row, or column)")
	names(cut.points) <- levels(group)

	if (is.null(layout))
		layout <- c(1, nlevels(group))
	obj <- histogram(~ log2(val) | group, data = df,
		layout = layout, n=64, as.table=TRUE,
		main = paste(levels(df$dname), collapse=" + "),
		panel = function(x, ..., subscripts)
		{
			panel.histogram(x,  ...)
			idx <- unique(group[subscripts])
			panel.abline(v=log2(cut.points[idx]), col=2)
		},
		scales=list(y=list(relation="free", rot=0)), ...
	)
	plot(obj)
	invisible(obj)
}
