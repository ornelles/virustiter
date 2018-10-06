#########################################################################################
# plotPlate
#
# display entire plate with lattice graphics (expects well values as "a03", etc.
# otherwise treats it as a 
#
#########################################################################################
#' Plot Results of an Entire Plate by Well or File
#'
#' Display an entire plate with lattice graphics grouped by well or file.
#'
#' @param df Annotated \code{data.frame} with imaging results.
#' @param size Number of wells for plate, if \code{NULL}, the size will be guessed.
#' @param cex Numeric value as factor by which to expand the plot symbol.
#' @param alpha Transparency factor applied to plotting symbols.
#' @param main Optional character string to serve as plot title.
#' @param invert.y A \code{logical} value to invert y coordinates.
#' @param ... Additional arguments handed to \code{xyplot} including subset.
#'
#' @import lattice
#'
#' @return
#'
#' The plot is generated and returned as an invisible \code{lattice} object.
#'
#' @export
#'
plotPlate <- function(df, size = NULL, cex = 1/2, alpha = 1/2, main = NULL,
	invert.y = TRUE, ...)
{
	require(lattice)

# default to well first
	if ("well" %in% names(df)) {
		byWell <- TRUE
		if (is.null(size))
			n <- nlevels(df$well)
		else {
			if (size < 1 | size > 384)
				stop("'size' must be between 1 and 384")
			n <- as.integer(size)
		}
		if (n > 96)	{rows <- 16; columns <- 24}     # 384-well plate
		else if (n > 48) {rows <- 8; columns <- 12} # 96-well plate
		else if (n > 24) {rows <- 6; columns <- 8}  # 48-well plate
		else if (n > 12) {rows <- 4; columns <- 6}  # 24-well plate
		else if (n > 6) {rows <- 3; columns <- 4}   # 12-well plate
		else if (n == 6) {rows <- 2; columns <- 3}  # 6-well plate
		else {rows <- 1; columns <- n}              # fewer than 6

	# create "fully populated" factor for well
		ww <- expand.grid(seq_len(columns), LETTERS[seq_len(rows)])
		ww <- apply(ww, 1, function(x) sprintf("%s%d",
			as.character(x[2]), as.numeric(x[1])))
			
	# ensure that well names are harmonized with full range of levels
		df$well <- well.info(df$well)$well
		df$well <- factor(df$well, levels = well.info(ww)$well)
		skip <- !levels(df$well) %in% unique(as.character(df$well))
	}
# else process by file
	else if ("file" %in% names(df)) {
		byWell <- FALSE
		n <- nlevels(df$file)
		flevels <- abbreviate(levels(df$file), 18, method = "both")
		df$file <- factor(df$file, levels = flevels) 
		rows <- columns <- ceiling(sqrt(n))
	}
	else
		stop("requires variable 'well' or 'file' in ", deparse(substitute(df)))

	if (is.null(main))
		main <- Sys.Date()
	if (!"positive" %in% names(df))
		df$positive <- FALSE
	if (invert.y)
		df$ym <- -df$ym
	if (byWell == TRUE) {
		obj <- xyplot(ym ~ xm | well, data = df, groups = positive, cex = cex,
			alpha = alpha, as.table = TRUE, aspect = "iso", layout = c(columns,rows),
			skip = skip, xlab = "", ylab = "", scales = list(draw = FALSE),
			main = main, ...)
	}
	else {
		obj <- xyplot(ym ~ xm | file, data = df, groups = positive, cex = cex,
			alpha = alpha, as.table = TRUE, aspect = "iso", layout = c(columns,rows),
 			xlab = "", ylab = "", scales = list(draw = FALSE), main = main, ...)
	}
	plot(obj)
	invisible(obj)
}
