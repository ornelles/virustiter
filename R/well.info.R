#' Extract ID, Prefix, Well, Row and Column Label
#' 
#' Create a uniform representation for a multi-titer well with an optional 
#' prefix extracted as prefix.
#' 
#' @param w label (tag) for the well (coerced to a character)
#' @param format character string as \code{\link{sprintf}} format for the
#'   column, default value of \code{NULL} pads the column value with zeros
#'   as needed such as \code{"\%02d"} for more than 9 columns
#' @param upper use upper case for row labels if \code{TRUE}
#' @param drop.levels drop unused levels in \code{row} and \code{column}
#'  otherwise the respective levels for \code{row} and \code{column} are
#'  the letters A to P and numbers 1 to 384.
#' 
#' @details 
#' 
#' This function will parse the \code{character} string in \code{w}
#' into an optional "harmonized" tag with an arbitrary character prefix,
#' (typically used to designate plate number), followed by a \strong{single}
#' character for the well ending with the column number formatted as per
#' the \code{sprintf} formatted \code{format}. The code will fail if the row 
#' value is not one of the 16 letters from A to P or if the column number
#' is not in the range of 1 to 384.
#'
#' The prefix must contain only alphanumeric characters (0-9, a-z, A-Z), the
#' underscore character or the hyphen.
#'
#' This function is also used to \emph{harmonize} well information generated
#' by different processes that may differ by the lower or upper case or the
#' use of leading zeros for the column number. 
#' 
#' @return 
#' 
#' A named list of length three or five. A value for \code{tag} and
#' \code{prefix} will be returned only if a prefix is identified. 
#' \itemize{
#'  \item tag, harmonized character representation of the entire well name
#'    including any prefix if present, otherwise this is not returned
#'  \item prefix, prefix as a character string if present, otherwise
#'    this is not returned
#'  \item well, harmonized well name as a factor
#'  \item row, row as a factor with levels A-P (or a-p)
#'  \item column, column number as a factor with levels 1-24
#' }
#' 
#' @examples
#'   well.info("a6") # default settings
#'   well.info("a6", format = "%d", upper = TRUE)
#'   well.info("Plate_2b92", format = "%05d")
#'   well.info("Plate_2b92", format = "%05d", drop.levels = FALSE)
#' 
#' @export
#' 
well.info <- function(w, format = NULL, upper = TRUE, drop.levels = TRUE)
{
# coerce to character, ensure logical values
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
  row <- substring(ww, 1, 1)  # row must be first position
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
