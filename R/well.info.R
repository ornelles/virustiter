#' Extract ID, Plate, Well, Row and Column Label
#' 
#' Create a uniform representation for a multi-titer well with an optional 
#' prefix extracted as plate 
#' 
#' @param w Label for the well (coerced to a character).
#' @param format Character string as \code{\link{sprintf}}
#'   format for the column (default of \code{"%02d"}). 
#' @param upper \code{logical} value to use upper case when \code{TRUE}.
#' 
#' @details 
#' 
#' This function will parse the \code{character} string in \code{w}
#' into a "harmonized" label with an optional arbitrary prefix used to
#' identify the plate, followed by a single character for the well,
#' followed by the  column number formatted as per \code{format}.
#' The code will fail if the row  value is not one of the 16 letters from
#' A to P or if the column number is not in the range of 1 to 384.
#'
#' The prefix must contain only alphanumeric characters (0-9, a-z, A-Z), the
#' underscore character, or the hyphen.
#'
#' This function is also used to \emph{harmonize} well information generated
#' by different processes that may differ by the lower or upper case or the
#' use of leading zeros for the column number. 
#' 
#' @return 
#' 
#' A named list of length four or five. 
#' \itemize{
#'  \item label, harmonized character representation of the entire label
#'    including any prefix
#' 	\item plate, prefix as a character string if present, otherwise
#'    this is not returned
#' 	\item well, harmonized well name as a factor
#' 	\item row, row as a factor with levels A-P (or a-p)
#' 	\item column, column number as a factor with levels 1-24
#' }
#' 
#' @examples
#'   well.info("a6") # default settings
#'   well.info("a6", format = "%d", upper = TRUE)
#'   well.info("Plate_2b92", format = "%05d")
#' 
#' @export
#' 
well.info <- function(w, format = "%02d", upper = TRUE)
{
# coerce to character
	w <- as.character(w)

# grep pattern for 'well' at end of string
	wellpat <- "[[:alpha:]][[:digit:]]+$"

# extract an optional prefix as everything but last letter followed by digits
	prefix <- sapply(strsplit(w, wellpat), "[", 1)

# extract information in well position: row and column
	ww <- sub(paste0("^.*(", wellpat, ")"), "\\1", w)
	ww <- toupper(as.character(ww))
	row <- substr(ww, 1, 1)	# row must be first position
	column <- substr(ww, 2, 12) # column must be 2nd position to end
	column <- as.integer(column)

# error checking
	if (any(grepl("[^[:alnum:]_-]", prefix)))
		stop("prefix must contain only letters, numbers, underscore or hyphen")
	if (!all(row %in% LETTERS[1:16]))
		stop("row must be a single character in the letters A - P  or a - p")
	if (any(column < 1 | column > 384))
		stop("column must be an integer between 1 and 384")

# format plate
	plate <- prefix

# format row
	if (upper == FALSE) { # convert back to lower case
		row <- tolower(row)
		row <- factor(row, levels = letters[1:16])
	}
	else
		row <- factor(row, levels = LETTERS[1:16])
	
# create harmonized well
	well <- paste(row, sprintf(format, column), sep = "")
	well <- factor(well)

# format column
	column <- factor(column, levels = 1:24)

# create harmonized label 
	label <- paste(plate, well, sep = "")

# assemble list with or without prefix value
	if (all(plate == ""))
		ans <- list(label = label, well = well, row = row, column = column)
	else
		ans <- list(label = label, plate = plate, well = well, row = row, column = column)
	return(ans)
}
