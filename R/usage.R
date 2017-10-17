#########################################################################################
# usage - display functional information
#
usage <- function() {
	txt <- c(
	" Usage:",
	"   df <- parseImages() # analyze images with EBImage",
	"      or",
	"   df <- readIJResults()   # read data from Fluorescent Cell Count v6 (ImageJ)",
	"",
	"   df  <- mergePdata(pd, df)  # merge with phenotype data in 'pd'",
	"   cut <- getCut(df)     # determine cutoff by control (or well, row, or column)",
	"   df  <- score(df, cut) # assign positive values from cutoff",
	"   res <- tally(df)      # tally positives and negatives and return data.frame",
	"   fm  <- getFit(obj)    # get model fit(s) from scored data (df) or from 'res'",
	"   cf  <- getTiter(fm)   # get value in units required for MOI of 1 and 95% CI",
	"", 
	" Support:",
	"   plotCut(df)    # calculate and show cutoff values with densityplot ",
	"   plotPlate(df)  # plot plate showing positives",
	"   plotWell(df, well) # plot each file in a well showing positives and sizes",
	"   plotHist(df)   # histogram of each well with optional cutoff values",
	"   plotFit(fm)    # plot fit(s) with calculated values using base graphics",
	"   plotOneFit(fm) # plot fit with options to adjust colors",
	"   addOneFit(fm)  # add best-fit line to existing base graph",
	"   getAIC(df, cut, by)  #evaluate fitted model(s) from df at cut values",
	"	displayPairs(f, dna = TRUE) # display image pairs in directory with 'f'",
	"",
	" Wrapper to automatically process results from ImageJ 'Results.txt' file",
	"   fitAndPlot(res, by)")
	ft <- tempfile()
	writeLines(txt, con = ft, sep = "\n")
	file.show(ft)
}
