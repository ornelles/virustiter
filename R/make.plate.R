#' Create data.frame with Tag, Prefix, Well, Row and Column
#' 
#' Create a data.frame representing a multi-well plate with an optional 
#' prefix.
#' 
#' @param n the number of wells (6, 12, 24, 48, 96, or 384)
#' @param prefix optional \code{character} vector of prefix values
#' @param digits optional \code{integer} value to zero-pad the column labels
#' @param upper \code{logical} flag to use upper case for row labels
#' @param byrow \code{logical} value specifying organization 
#' 
#' @details 
#' 
#' Generate a data.frame suitable for typical phenodata. The \code{tag} 
#' is formed by replicating the prefix and pasting it in
#' front of each well value.
#'
#' @return 
#' 
#' A data.frame of length three or five with the following items. Note that 
#' values for \code{tag} and \code{prefix} will be returned only if
#' a non-NULL value for prefix is provided. 
#' \itemize{
#'  \item tag, harmonized character representation of the entire well name
#'    including any prefix, if a prefix was provided
#'  \item prefix, prefix as a character string, if provided
#'  \item well, harmonized well name as a factor
#'  \item row, row as a factor 
#'  \item column, column number as a factor
#' }
#' 
#' @examples
#'   make.plate(6)
#'   make.plate(6, upper = FALSE)
#'   make.plate(12, 1:3)
#'   make.plate(24, 1:6)$tag # typical use to generate "OTE" labels
#' 
#' @export
#' 
make.plate <- function(n, prefix = NULL, digits = NULL, upper = TRUE,
  byrow = TRUE)
{
# process arguments
  n <- as.integer(n)
	stopifnot(length(n) == 1)
  if (!is.null(prefix)) prefix <- as.character(prefix)
  if (!is.null(digits)) digits <- as.integer(digits)[1]
  upper <- as.logical(upper)
  byrow <- as.logical(byrow)

# assign number of rows (and columns) from desired plate size
  if (n == 6) nr <- 2
  else if (n == 12) nr <- 3
  else if (n == 24) nr <- 4
  else if (n == 48) nr <- 6
  else if (n == 96) nr <- 8
  else if (n == 384) nr <- 16
  else stop("'n' must be 6, 12, 24, 48, 96 or 384")
  nc <- n/nr

# determine padding needed for column number
  if (is.null(digits)) digits <- ceiling(log10(nc))

# determine case for row
  rowvals <- if (upper == TRUE) LETTERS[1:16] else letters[1:16]
    
# create format string for well
  fmt <- paste0("%s%0", digits, "d")
  if (byrow == TRUE)
    v <- expand.grid(column = 1:nc, row = rowvals[1:nr])[,2:1]
  else
    v <- expand.grid(row = rowvals[1:nr], column = 1:nc)
  well <- apply(v, 1, function(x) sprintf(fmt, x[1], as.numeric(x[2])))
    
# assemble well, row and column
  well <- factor(well)
  column <- factor(v$column, levels = 1:nc)
  row <- factor(v$row, levels = rowvals[1:nr])
  ans <- data.frame(well = well, row = row, column = column) 

# process optional prefix and assemble data.frame
  if (!is.null(prefix)) {
    prefix <- rep(prefix, each = length(well))
    tag <- paste0(prefix, as.character(well))
    ans <- cbind(tag = tag, prefix = prefix, ans)
  }
# return data.frame 
  return(ans)
}
