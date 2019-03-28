#' Determine Background Value to Define Positive Samples
#'
#' Determine the optimum background value between "positive" and "negative"
#' values in the argument \code{param}. The group used to determine the
#' background is specified by the argument \code{by}.  
#'
#' @param df Annotated \code{data.frame} with fluorescent values to evaluate.
#' @param by Character string identifying the grouping factor in \code{'df'}
#'   used to determine background values. Values of "row" or "column" will
#'   split the data by row or column before identifying background values.
#'   If missing or if specified as "row" or "column", values associated with
#'   \code{type == "control"} or \code{moi == 0} or \code{x == 0} will
#'   be used to determine the background.
#' @param param Variable name in \code{df} as a character string to evaluate, 
#'   typically \code{"mfi"} or \code{"y"}.
#' @param mult Muliplier constant passed to \code{\link{findBgnd}}.
#' @param log \code{logical} flag passed to \code{\link{findBgnd}} to use
#'   log-transformed values.
#'
#' @details
#' #'
#' The value between positive and negative values in \code{param} will be 
#' determined according to \code{by}. If this value is "control" or is missing,
#' all values identified as \code{type == "control"} or with \code{x/moi == 0}
#' will be treated as background. If this value is "row" or "column", 
#' any values with \code{type == "control"} or with \code{x/moi == 0} in each
#' row or column will be used to define the background for that row or column.
#' If zero moi values are not present, the background will be determined by
#' Otsu's method for the groups defined by the character string in \code{by}.
#' Typically this would be \code{"well"} or \code{"file"} but can be any factor 
#' variable in the data.frame \code{df}. The background value will be determined
#' by the logic in \code{\link{findBgnd}}.
#' 
#' The annotated data frame must have the variable identified in \code{param}
#' and, if 'by' is missing, a variable named either \code{"x"} or \code{"moi"}.
#' If 'by' is provided, this must exist as a factor in the in the argument
#' \code{df}.
#'
#' @return
#'
#' A named numeric vector of background values.
#'
#' @examples
#' # Subset of data by_stack
#'   f <- system.file("extdata", "by_stack/file005.tif", package = "virustiter")
#'   img <- getImages(f)
#'   v <- parseImages(img)
#'
#' # Get background value for grouping value of "file"
#'   getBgnd(v, by = "file")
#'
#' @import EBImage
#'
#' @export
#'
getBgnd <- function(df, by, param = "mfi", mult = 2.5, log = TRUE)
{
	if (missing(df)) {
		usage <- c("getBgnd examples:",
			'  getBgnd(df, by = "control", param = "mfi", mult = 2.5, log = TRUE)',
			'  getBgnd(df) # same as above',
			'  getBgnd(df, "row", mult = 3, log = FALSE)')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}
	if (!is.data.frame(df))
		stop(deparse(substitute(df)), " must be a data.frame")

# parse 'by' argument
	if (missing(by) || by == "control")
		by <- "control"
	else if(!is.character(by))
		stop ("'", deparse(substitute(by)), "' must be a character string")
	else {
		index <- pmatch(by, names(df))
		if (is.na(index))
			stop('"', by, '" is not in data frame')
		else
			by <- names(df)[index]
	}

# parse 'param' argument
	if (!param %in% names(df))
		stop("'", deparse(substitute(param)), "' is not in data frame")

# assign name to 'moi'
	if ("moi" %in% names(df)) moi <- "moi"
	else if ("x" %in% names(df)) moi <- "x"
	else moi <- ""
		
# extract the appropriate grouping variable and mfi data in 'temp'
	if (by == "control") {
		if (any(df$type == "control"))# priority over zero moi values
			temp <- data.frame(g = TRUE, y = df[[param]][df$type == "control"])
		else if (any(df[[moi]] == 0))
			temp <- data.frame(g = TRUE, y = df[[param]][df[[moi]] == 0])
		else
			stop("\n", "If 'by' is missing or is specified as \"control\"",
				"\n", "type == \"control\" or 'x/moi' == 0 values must be present")
	}
	else if (by %in% c("column", "row")) { # try to select moi == 0 cases
		spl <- split(df, df[[by]])
		idx <- lapply(spl, function(v) which(v[[moi]] == 0))
		sel <- which(lengths(idx) == 0)
		idx[sel] <- TRUE
		temp <- do.call(rbind, Map(function(v, s) v[s, ], spl, idx))[c(by, param)]
		names(temp) <- c("g", "y")
		levels(temp$g) <- factor(temp$g, levels = unique(temp$g))
	}
	else {
		temp <- df[c(by, param)]
		names(temp) <- c("g", "y")
		levels(temp$g) <- levels(df[[by]])
	}
	res <- aggregate(y ~ g, temp, function(v) findBgnd(v, mult = mult, log = log))
	ret <- c(res[["y"]])
	names(ret) <- if(by == "control") rep("control", length(ret)) else res$g
	return(ret)
}
