#########################################################################################
# plotHist
#
# display histogram on well by well basis with lattice graphics, display background
# 'cut' values be a single value or should be a named vector from getCut()
#
#########################################################################################

plotHist <- function(df, cut=NULL, layout=NULL, ...) {
	if (missing(df)) {
		usage <- c("plotHist examples:",
			'  plotHist(df)      ## uses "positive" values in df',
			'  plotHist(df, cut) ## where cut was from getCut()')
		cat(usage, sep="\n")
		return(invisible(NULL))
	}
	library(lattice)
	if (is.null(cut))					# use existing 'positive' assignment
		cut.points <- with(df, tapply(val, list(positive, well), max))[1,]
	else if (length(cut) == 1)			# single value
		cut.points <- rep(cut, nlevels(df$well))
	else if (all(names(cut) %in% levels(df$well)))
		cut.points <- cut
	else if (all(names(cut) %in% levels(df$row)))
		cut.points <- rep(cut, each=nlevels(df$column))
	else if (all(names(cut) %in% levels(df$column)))
		cut.points <- rep(cut, nlevels(df$row))
	else
		stop("cut must be null, a single number, or a named vector (well, row, or column)")
	names(cut.points) <- levels(df$well)

	if (is.null(layout))
		layout <- c(1, nlevels(df$well))
	obj <- histogram(~ log2(val) | well, data = df,
		layout = layout, n=64, as.table=TRUE,
		main = paste(levels(df$dname), collapse=" + "),
		panel = function(x, ..., subscripts)
		{
			panel.histogram(x,  ...)
			well <- unique(df$well[subscripts])
			panel.abline(v=log2(cut.points[well]), col=2)
		},
		scales=list(y=list(relation="free", rot=0)), ...
	)
	plot(obj)
	invisible(obj)
}
