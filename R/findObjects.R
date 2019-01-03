#' Find Objects - NOT QUITE
#' 
#' Identify objects in data.frame produced by \code{parseImages()} 
#' as specified by the expression \code{expr}. 
#' 
#' @param expr An \code{expression} evaluated in \code{'data'} that returns
#'   a logical vector identifying objects (rows) in \code{'data'}.
#' @param mask An optional \code{Image} object in \code{Grayscale} mode or
#'   an array containing object masks. Object masks are sets of pixels with
#'   the same unique integer value.
#' @param data A \code{data.frame} produced by \code{parseImages()} in which
#'   \code{'expr'} will be evaluated.
#' @param invert A \code{logical} value to invert the selection determined
#'   by \code{'expr'}.
#' 
#' @return
#'
#' An integer \code{Image} mask or list of integer \code{Image} masks with the
#' objects selected by \code{'expr'} if \code{mask} was provided otherwise
#' a nested list of integers identifying the selected objects.
#' 
#' @examples
#'   x <- getImages(system.file("extdata", "by_folder/b4", package = "virustiter"))
#'   nm0 <- nucMask(x$nuc)
#'   df0 <- parseImages(x$nuc, x$tgt, nm0)
#'   nm1 <- findObjects(area > 190 & area < 220, df0, nm0)
#'   a <- Map(paintObjects, nm1, lapply(x$nuc, function(v) toRGB(normalize(v))))
#'   plot(a$b4, all = TRUE)
#'
#' @import EBImage
#'
#' @export
#'
findObjects <- function(expr, data, mask, invert = FALSE)
{
	if (missing(expr) || missing(data)) {
		usage <- c("findObjects examples:",
			'  findObjects(area > 200, df) # identify objects with area > 200',
			'  findObjects(positive == TRUE, df, mask) # mark positive objects in mask')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}
# check 'data' argument
	if (!is(data, "data.frame"))
		stop("'", deparse(substitute(data)), "' must be a data.frame")
	dnames <- names(data)
	if (!("frame" %in% dnames & ("well" %in% dnames | "file" %in% dnames)))
		stop('"frame" and either "well" or "file" must be in ', deparse(substitute(data)))

# evaluate expression and assign to original data.frame
	sel <- eval(substitute(expr), data)
	vname <- tail(make.unique(c(names(data), "var")), 1)
	data[[vname]] <- sel

# split 'var' according to 'well/file' and 'frame' in 'data'
	if ("well" %in% names(data))
   group <- "well"
	else if ("file" %in% names(data))
		group <- "file"
	spl.1 <- split(data, data[[group]], drop = FALSE)
	spl.2 <- lapply(spl.1, function(v) split(v[[vname]], v[["frame"]], drop = TRUE))
	neg <- lapply(spl.2, function(v) lapply(v, function(x) which(x == FALSE)))
	pos <- lapply(spl.2, function(v) lapply(v, function(x) which(x == TRUE)))

# process mask
	if (missing(mask))
		ans <- if(invert) neg else pos
	else if (is(mask, "Image")) {
		pos <- unlist(pos, recursive = FALSE)
		neg <- unlist(neg, recursive = FALSE)
		ans <- if(invert) rmObjects(mask, pos) else rmObjects(mask, neg)
	}
	else if (is(mask, "list"))
		ans <- if(invert) Map(rmObjects, mask, pos) else Map(rmObjects, mask, neg)
	else
		stop("unable to use combination of 'data' and 'mask'")
	return(ans)
}
