#' Extract ID, Prefix, Well, Row and Column Label
#' 
#' Create a uniform representation for a multi-titer well with an optional 
#' prefix extracted as prefix \emph{or} generate well labels.
#' 
#' @param w label (tag) for the well (coerced to a character) \emph{or} a
#'  single numeric value of 6, 12, 24, 48, 96, or 384 to generate well labels.
#' @param format Character string as \code{\link{sprintf}} format for the
#'   column, default value of \code{NULL} pads the column value with zeros
#'   as needed such as \code{"\%02d"} for more than 9 columns. 
#' @param upper to use upper case when \code{TRUE}.
#' @param drop.levels \code{logical} value to drop unused levels in
#'   \code{row} and \code{column}.
#' @param byrow \code{logical} value when generating well labels.
#' 
#' @details 
#' 
#' This function will parse the \code{character} string in \code{w}
#' into an optional "harmonized" tag with an arbitrary prefix used to
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
#' If \code{w} is an integer, well labels ("A1", "A01", etc) will be generated
#' for the requested number of wells (6, 12, 24, 48, 96 or 384) in the direction
#' indicated by \code{byrow}.
#'
#' @return 
#' 
#' A character vector of well names or a named list of length three or five.
#' Values for \code{tag} and
#' \code{prefix} will be returned only if a plate prefix is used. 
#' \itemize{
#'  \item tag, harmonized character representation of the entire well name
#'    including any prefix if present, otherwise this is not returned
#' 	\item prefix, prefix as a character string if present, otherwise
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
#'   well.info("Plate_2b92", format = "%05d", drop.levels = FALSE)
#'   well.info(12)
#' 
#' @export
#' 
well.info <- function(w, format = NULL, upper = TRUE, drop.levels = TRUE,
	byrow = TRUE)
{
# check for special case of w as numeric
	byrow <- as.logical(byrow)
	if (is(w, "numeric") && length(w) == 1) {
		if (w == 6) nr <- 2
		else if (w == 12) nr <- 3
		else if (w == 24) nr <- 4
		else if (w == 48) nr <- 6
		else if (w == 96) nr <- 8
		else if (w == 384) nr <- 16
		else
			stop("'w' must be 6, 12, 24, 48, 96 or 384 to generate well labels")
		nc <- w/nr
		fmt <- if (w > 48) "%s%02d" else "%s%d" 
		if (byrow == TRUE) v <- expand.grid(column = 1:nc, row = LETTERS[1:nr])[,2:1]
		else v <- expand.grid(row = LETTERS[1:nr], column = 1:nc)
		return(apply(v, 1, function(x) sprintf(fmt, x[1], as.numeric(x[2]))))
	}
		
# otherwise coerce to character, ensure logical values
	w <- as.character(w)
	upper <- as.logical(upper)
	drop.levels <- as.logical(drop.levels)

# grep pattern for 'well' at end of string
	wellpat <- "[[:alpha:]][[:digit:]]+$"

# extract an optional prefix as everything but last letter followed by digits
	prefix <- sapply(strsplit(w, wellpat), "[", 1)

# extract information in well position: row and column
	ww <- sub(paste0("^.*(", wellpat, ")"), "\\1", w)
	ww <- toupper(as.character(ww))
	row <- substring(ww, 1, 1)	# row must be first position
	column <- substring(ww, 2) # column must be 2nd position to end
	column <- as.integer(column)

# assign format string
	if (is.null(format)) {
		digits <- ceiling(log10(max(column)))
		format <- paste0("%0", digits, "d")
	}

# error checking
	if (any(grepl("[^[:alnum:]_-]", prefix)))
		stop("prefix must contain only letters, numbers, underscore or hyphen")
	if (!all(row %in% LETTERS[1:16]))
		stop("row must be a single character in the letters A - P  or a - p")
	if (any(column < 1 | column > 384))
		stop("column must be an integer between 1 and 384")

# format row
	if (upper == FALSE) { # convert back to lower case
		row <- tolower(row)
		row <- factor(row, levels = letters[1:16])
	}
	else
		row <- factor(row, levels = LETTERS[1:16])

# create harmonized well as factor with default levels
	well <- paste(row, sprintf(format, column), sep = "")
	well <- factor(well)

# format column
	column <- factor(column, levels = 1:384)

# create harmonized tag 
	tag <- paste(prefix, well, sep = "")

# drop levels from row and column?
	if (drop.levels) {
		row <- droplevels(row)
		column <- droplevels(column)
	}

# assemble list with or without prefix value
	if (all(prefix == ""))
		ans <- list(well=well, row=row, column=column)
	else
		ans <- list(tag=tag, prefix=prefix, well=well, row=row, column=column)
	return(ans)
}
