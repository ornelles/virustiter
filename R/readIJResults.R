#########################################################################################
# readIJResults
#
# read data from ImageJ results file and extract infection parameters into a data.frame
# This program expects output from the ImageJ code "Multiwell Fluorescent Cell Count
# Expected values from ImageJ are: Area, Mean gray value, Center of Mass, Display label,
# and Inverted Y coordinates. Label must be of the form "xnn/filename" where xnn is the well
# such as "A01" or "b2". The ImageJ results file should have the multiplicity as an x value
# and a unit indicator (unit) coded as 1, 2, 3 or 4 for
# "VP per cell", "IU per cell", "ul per cell", "ml per cell"
#
# Arguments
#	f			optional ImageJ results text file
#	SEP			character used to separate directory and filename in ImageJ results file
#
# Dependencies
#	well.info()	parses well information
#
# Returns invisible raw data with extracted well, row, column, and positive (as FALSE)
#
#########################################################################################

readIJResults <- function(f, SEP = NULL)
{
# coded vector of units from ImageJ, note that the "per well" units were
# converted to per cell in the ImageJ code

	unitList <- c("VP", "IU", "ul", "ml")	# used in "1 IU = 2.5 VP, etc.

# read data from ImageJ results file and adjust ImageJ names
	if (missing(f))
		f <- file.choose()
	if (is.null(SEP))
		SEP <- .Platform$file.sep
	df <- read.table(f, header=TRUE)
	names(df) <- tolower(names(df))
	sel <- which(names(df)=="mean")
	names(df)[sel] <- "val"
	sel <- which(names(df)=="label")
	names(df)[sel] <- "file"


# extract unit string from first entry in Image J Results file
	if("unit" %in% names(df))
		unit <- unitList[df$unit[1]]
	else
		unit <- "none"

# extract well, row, column and file name information such as "a01/file001.tif"
	txt <- strsplit(as.character(df$file), SEP)
	if (!all(sapply(txt, length) == 2))
		stop("encountered file label with more than one separator, ", SEP)
	txt <- unlist(txt)
	v <- txt[seq(1,length(txt),2)]
	row <- well.info(v)$row				# row, checks for only letters
	column <- well.info(v)$column		# column, checks for sensible numbers
	well <- well.info(v)$well			# well, a01 format returned
	dname <- basename(dirname(f))		# parent directory name
	fname <- txt[seq(2,length(txt),2)]	# Image file name
	fname <- gsub("\\..*$", "", fname)	# file name without extension

# assemble first part of data.frame
	df <- data.frame(df[1], fname=fname, directory=dname, well=well, row=row,
			column=factor(column), df[-1], unit)

# Extract x values from Image J results file and verify that each well is unique.
	if (!("x" %in% names(df)))
		stop("The ImageJ results file must have 'x' values")
	xx <- with(df, tapply(x, well, unique))		# xx holds an array of x values
	tally <- sapply(xx, length)					# count unique 'x' values
	if (length(tally) == 1)
		warning("only a single x value was found: ", xx[1])
	if(!all(tally == 1))
		warning("multiple x values found: ", names(which(tally > 1)))

# assign type to each well as standard or control if x==0
	type <- rep("standard", nrow(df))
	control.wells <- names(which(xx==0))
	type[df$well %in% control.wells] <- "control"

# assemble and return data.frame, change "mean" to "val"
	df <- data.frame(df, type=type)
	return(df)
}
