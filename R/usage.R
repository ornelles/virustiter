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
  "   img <- getImages()     # read paired images with EBImage",
  "   df <- parseImages(img) # extract nuclear information and target mfi",
	"",
	"   pd <- data.frame(well = levels(df$well), moi = moi, unit = unit) ...or...",
	"   pd <- data.frame(file = levels(df$file), moi = moi, unit = unit)",
	"",
	"   df  <- mergePdata(pd, df) # merge with phenotype data in 'pd'",
	"   bg <- getBgnd(df)     # determine background by control (or well, row, or column)",
	"   df  <- score(df, bg)  # assign positive values from background",
	"   res <- tally(df)      # tally positives and negatives and return data.frame",
	"   fm  <- getFit(res)    # get model fit(s) from either res or scored data frame",
	"   cf  <- getTiter(fm)   # get titer infectious units per ml +/- 95% CI",
	"", 
	" Support:",
	"   checkImages(path)   # check (and display) paired images", 
	"   plotHist(df)        # histogram of each well with optional background values",
	"   plotDens(df)        # calculate and show background values with densityplot ",
	"   plotPlate(df)       # plot plate showing positives",
	"   plotWell(df, well)  # plot each cell in a well showing positives and sizes",
	"   plotFit(fm)         # plot fit(s) with calculated values using base graphics",
	"   plotOneFit(fm)      # plot fit with options to adjust colors",
	"   addOneFit(fm)       # add best-fit line to existing base graph",
	"   getAIC(df, bg)      # evaluate fitted model(s) from df at bg values",
	"   getTiter(fm)        # get titer as IU per ml (if volume units were used)",
  "   getEC63(fm)         # get volume required for 1 infectious units +/- 95% CI",
	"   nucMask(dapi)       # extract nuclear mask from dapi image(s) or file(s)",
	"   trimMask(mask)      # remove objects based on size from mask",
	"   cellMask(mask)      # expand a nuclear mask into a cell mask",
  "   edgeObjects(mask)   # identify objects near the edge of a mask",
	"   getVal(mask, tgt)   # extract single 'computeFeatures' value (b.mean) from tgt",
	"   p2p()               # interactively measure point-to-point distances",
	"")
	ft <- tempfile()
	writeLines(txt, con = ft, sep = "\n")
	file.show(ft, title = "Information on 'virustiter' package")
}
