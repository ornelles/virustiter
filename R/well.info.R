#' Extract Prefix, Row and Column from Well Label
#' 
#' Create a uniform representation for a well with an optional 
#' prefix. 
#' 
#' @param w Name of the well (coerced as a character).
#' @param format Character string as \code{\link{sprintf}} format for the column portion.
#' @param upper \code{logical} value to use upper case when \code{TRUE}.
#' 
#' @details 
#' 
#' Parse the value in \code{w} into an optional arbitrary prefix followed by a 
#' single character for the well, followed by the column number formatted
#' as per \code{format}. The code will fail if the row value is not one of the 26
#' letters or if the column number is not in the range of 1 to 384. 
#' 
#' @return 
#' 
#' A named list of length three or four:
#' \itemize{
#' 	\item prefix, the prefix if present; otherwise this is not returned
#' 	\item well, uniform character representation of the well
#' 	\item row, single character for the row
#' 	\item column, column number as a factor
#' }
#' 
#' @examples
#'   well.info("2a6") # default settings
#'   well.info("2a6", format = "%d", upper = TRUE)
#'   well.info("b92", format = "%05d")
#' 
#' @export
#' 
well.info <- function(w, format = "%02d", upper = TRUE)
{
# coerce to character
	w <- as.character(w)

# extract optional prefix as everything but last letter followed by digits
	prefix <- sapply(strsplit(w, "[[:alpha:]][[:digit:]]+$"), "[", 1)

# extract well, row and column
	ww <- sub(".*([[:alpha:]][[:digit:]]+$)", "\\1", w)
	ww <- toupper(as.character(ww))
	row <- substr(ww, 1, 1)	# row must be first position
	column <- substr(ww, 2, 12) # column must be 2nd position to end
	column <- as.integer(column)

# error checking
	if (!all(row %in% LETTERS[1:26]))
		stop("bad row value")
	if (any(column < 1 | column > 384))
		stop("bad column value")
	if (upper == FALSE)
		row <- tolower(row)
	well <- paste(row, sprintf(format, column), sep = "")
	column <- factor(column, levels = 1:max(column))
	if (all(prefix == ""))
		ans <- list(well = well, row = row, column = column)
	else
		ans <- list(prefix = prefix, well = well, row = row, column = column)
	return(ans)
}
