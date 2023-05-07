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
	"   pd <- data.frame(file = levels(df$file), moi = moi, unit = unit)",
	"      Or for a single plate organized by well:",
	"   pd <- data.frame(well.info(unique(df$well)), moi = moi, unit = unit)",
	"      Or if multiple plates are used:",
	"   pd <- data.frame(well.info(unique(df$labels)), moi = moi, unit = unit)",
	"",
	"   df  <- mergePdata(pd, df) # merge with phenotype data in 'pd'",
	"   bg  <- getBgnd(df)   # determine background by control (or well, row, or column)",
	"   df  <- score(df, bg) # assign positive values from background",
	"   res <- tally(df)     # tally positives and negatives and return data.frame",
	"   fm  <- getFit(res)   # get model fit(s) from either res or scored data frame",
	"   cf  <- getTiter(fm)  # get titer infectious units per ml +/- 95% CI",
	"", 
	" Support:",
	"   checkImages(path)   # check (and display) paired images", 
  "   list.images(path)   # list image files in the given path",
	"   getZero(x)          # find most common non-zero pixel in each image frame",
	"   setZero(x)          # rescale each image frame to the most common zero value",
	"   nucMask(dapi)       # extract nuclear mask from dapi image(s) or file(s)",
	"   trimMask(mask)      # remove objects based on size from mask",
	"   cellMask(mask)      # expand a nuclear mask into a cell mask",
  "   edgeObjects(mask)   # identify objects near the edge of a mask",
  "   findObjects(expr, df) # find objects in data.frame identified by expr",
	"   getVal(mask, ref)   # extract one 'computeFeatures' value using mask and ref",
	"   p2p()               # interactively measure point-to-point distances",
	"   plotHist(df)        # histogram of each well with optional background values",
	"   plotDens(df)        # calculate and show background values with densityplot",
	"   plotPlate(df)       # plot plate showing positives",
	"   plotWell(df, well)  # plot each cell in a well showing positives and sizes",
	"   plotFit(fm)         # plot fit(s) with calculated values using base graphics",
	"   addOneFit(fm)       # add best-fit line to existing base graph",
	"   getAIC(df, bg)      # evaluate fitted model(s) from df at bg values",
	"   getTiter(fm)        # get titer as IU per ml (if volume units were used)",
  "   getEC63(fm)         # get volume required for 1 infectious units +/- 95% CI",
	"   getShift(mask, tgt) # get optimal x-y shift to align tgt with (nuc) mask",
	"   translate(x, v)     # apply (optimal) x-y shift in v to image x",
	"")
	ft <- tempfile()
	writeLines(txt, con = ft, sep = "\n")
	file.show(ft, title = "Information on 'virustiter' package")
}
