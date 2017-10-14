#########################################################################################
#
# fitAndPlot
#
# wrapper to process from results object or ImageJ Results.text file, as a sequential
# (single plate) or by columns or rows. Control wells (moi == 0) are expected to be present.
# In addition to plotting and printing the titer, the adjusted data frame is returned
# as an invisible ojbect
#
#########################################################################################

fitAndPlot <- function(f, by=c("sequential", "column", "row"), plot.plate=FALSE)
{
	if (missing(f)) f <- file.choose()
	by <- match.arg(by)
	df <- readIJResults(f)
	if (by == "sequential") cut <- getCut(df, "control")
	else cut <- getCut(df, by)
	df <- score(df, cut)
	res <- tally(df)
	fm <- getFit(res, by)
	cf <- getTiter(fm)
	if (plot.plate) (plotPlate(df))
	n <- length(fm)
	if (n == 30) n <- 1 # single fit versus list of fits
	nr <- min(n, 4)
	nc <- ceiling(n/4)
	dev.new()
	par(mfrow=c(nr, nc))
	plotFit(fm)
	print(cf)
	invisible(df)
}

