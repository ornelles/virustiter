#' Show Usage of \code{virustiter} Package
#'
#' Display succinct summary of functions and work flow in \code{virustiter}.
#'
#' @return
#'
#' No value returned, a window is opened with helpful information.
#'
#' @export
#'
usage <- function() {
	txt <- c(
	" Workflow:",
	"   df <- parseImages()   # read paired images with EBImage ...or...",
	"   df <- readIJResults() # read data from Fluorescent Cell Count (ImageJ)",
	"",
	"   pd <- data.frame(well = levels(df$well), moi = moi, unit = unit) ...or...",
	"   pd <- data.frame(file = levels(df$file), moi = moi, unit = unit)",
	"",
	"   df  <- mergePdata(pd, df) # merge with phenotype data in 'pd'",
	"   cut <- getCut(df)     # determine cutoff by control (or well, row, or column)",
	"   df  <- score(df, cut) # assign positive values from cutoff",
	"   res <- tally(df)      # tally positives and negatives and return data.frame",
	"   fm  <- getFit(res)    # get model fit(s) from either res or scored data frame",
	"   cf  <- getTiter(fm)   # get value in units required for MOI of 1 and 95% CI",
	"", 
	" Support:",
	"   plotDens(df)        # calculate and show cutoff values with densityplot ",
	"   plotHist(df)        # histogram of each well with optional cutoff values",
	"   plotPlate(df)       # plot plate showing positives",
	"   plotWell(df, well)  # plot each cell in a well showing positives and sizes",
	"   plotFit(fm)         # plot fit(s) with calculated values using base graphics",
	"   plotOneFit(fm)      # plot fit with options to adjust colors",
	"   addOneFit(fm)       # add best-fit line to existing base graph",
	"   getAIC(df, cut)     # evaluate fitted model(s) from df at cut values",
	"   nucMask(dapi)       # extract nuclear mask from dapi image(s) or file(s)",
	"   trimMask(mask)      # remove objects based on size from mask",
	"   cellMask(mask)      # expand a nuclear mask into a cell mask",
	"   p2p()               # interactively measure point-to-point distances",
	"   pnpoly(p, v)        # test if points in p are within polygon (v)",
	"")
	ft <- tempfile()
	writeLines(txt, con = ft, sep = "\n")
	file.show(ft, title = "Information on 'virustiter' package")
}
