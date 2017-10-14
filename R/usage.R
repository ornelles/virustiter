#########################################################################################
# usage - display functional information
#
usage <- function() {
	txt <- c(
	" Usage:",
	"   df <- parseImages()   # read with EBImage",
	"      or",
	"   df <- readData()      # read data from Fluorescent Cell Count v6 (ImageJ)",
	"",
	"   df  <- mergePdata(pd, df)  # optional merge with phenoData in 'pd'",
	"   cut <- getCut(df)     # determine cutoff by control (or well, row, or column)",
	"   df  <- score(df, cut) # assign positive values from cutoff",
	"   res <- tally(df)      # tally positives and negatives and return data.frame",
	"   fm  <- getFit(res)    # get model fit(s)",
	"   cf  <- getTiter(fm)   # get reciprocal of titer and 95% confidence intervals",
	"", 
	" Support:",
	"   plotCut(df)   # calculate and show cutoff values with densityplot ",
	"   plotPlate(df) # plot plate showing positives",
	"   plotWell(df, well) # plot each file in a well showing positives and sizes",
	"   plotHist(df)  # histogram on well-by-well basis with optional cutoff values",
	"   plotFit(fm)   # plot fit(s) with calculated values using base graphics",
	"   plotOneFit(fm)  # plot fit with options to adjust colors",
	"   addOneFit(fm) # add best-fit line to existing base graph",
	"   getAIC(df, cut, by)  #evaluate fitted model(s) from df at cut values",
	"",
	" Wrapper to automatically process results data frmae or ImageJ 'Results.txt' file",
	"   fitAndPlot(res, by)")
	ft <- tempfile()
	writeLines(txt, con = ft, sep = "\n")
	file.show(ft)
}