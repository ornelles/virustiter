#' Determine Cutoff Value between Positive and Negative
#'
#' Determine the cutoff value between "positive" and "negative" values in
#' argument \code{param}. The background values are specified by \code{by} as
#' control values or by row, column or file.  
#'
#' @param df Annotated \code{data.frame} with fluorescent values to evaluate.
#' @param by Character string identifying the group in which to seek cutoff
#'   values \code{("control", "file", "well", "row" or "column")}.
#' @param param Variable name as character string in \code{df} to evaluate, 
#'   typically \code{"y"}.
#' @param mult Muliplier constant passed to \code{findBgnd()}.
#' @param log \code{logical} flag passed to \code{findBgnd()} to use
#'   log-transformed values.
#'
#' @details
#'
#' The cutoff value between positive and negative values in 
#' \code{param} will be determined according to \code{by}. If this
#' value is "control", all values identified as 
#' \code{type == "control"} or with \code{x == 0} will be considered to be 
#' background. Otherwise, a cutoff will be determined by Otsu's method
#' for the groups identified identified in \code{by} ("file", "well",
#' "column" or "row"). The cutoff value will be determined by the logic
#' in \code{findBgnd()}.
#' 
#' The annotated data frame must have the variable identified in 
#' \code{param} and a variable named \code{x} as well as an appropriate
#' grouping value: either \code{file, well, row} or \code{column}
#' as required. 
#'
#' @return
#'
#' A named numeric vector of cutoff values.
#'
#' @examples
#' # Subset of data by_stack
#'   f <- system.file("extdata", "by_stack/file005.tif", package = "virustiter")
#'   file.copy(f, tempdir(), overwrite = TRUE) # need files to reside in directory
#'   v <- parseImages(tempdir(), pattern = "file005")
#'
#' # Get cutoff for grouping value of "file"
#'   getCut(v, by = "file")
#'
#' @import EBImage
#'
#' @export
#'
getCut <- function(df, by = c("control", "file", "well", "row", "column"),
	param = "mfi", mult = 5, log = TRUE)
{
	if (missing(df)) {
		usage <- c("getCut examples:",
			'  getCut(df, by = "control", param = "mfi", mult = 5, log = TRUE)',
			'  getCut(df) # same as above',
			'  getCut(df, "row", mult = 3, log = FALSE)')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}
	if (!is.data.frame(df))
		stop(deparse(substitute(df)), " must be a data.frame")

# parse arguments and perform error checking
	by <- match.arg(by)
	if (!param %in% names(df))
		stop("'", param, "' not in data frame")

# create local copy of relevant data
	if (by == "control") { # special case of by "control"
		if (any(df$type == "control"))
			temp <- data.frame(g = TRUE, y = subset(df, type == "control")[[param]])
		else if (any(df$x == 0))
			temp <- data.frame(g = TRUE, y = subset(df, x == 0)[[param]])
		else
			stop('require type == "control" or x values of 0 to use by = "control" option')
	}
	else {
		temp <- df[c(by, param)]
		names(temp) <- c("g", "y")
	}
	res <- aggregate(y ~ g, temp, function(v) findBgnd(v, mult = mult, log = log))
	ret <- res$y
	names(ret) <- if(by == "control") "control" else res$g
	return(ret)
}
