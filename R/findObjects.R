#' Find Objects
#' 
#' Identify objects in data.frame produced by \code{parseImages()} 
#' as specified by the expression \code{expr}. 
#' 
#' @param expr An \code{expression} evaluated in \code{'df'} that returns
#'   a logical vector identifying objects (rows) in \code{'df'}.
#' @param mask An optional \code{Image} object in \code{Grayscale} mode or
#'   an array containing object masks. Object masks are sets of pixels with
#'   the same unique integer value.
#' @param df A \code{data.frame} produced by \code{parseImages()} in which
#'   \code{'expr'} will be evaluated.
#' @param invert A \code{logical} value to invert the selection determined
#'   by \code{'expr'}.
#' 
#' @return
#'
#' Either an integer \code{Image} mask or list of integer \code{Image} masks
#' with the objects selected by \code{'expr'} if \code{mask} was provided
#' otherwise a nested list of integers identifying the selected objects.
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
findObjects <- function(expr, df, mask = NULL, invert = FALSE)
{
	if (missing(expr) || missing(df)) {
		usage <- c("findObjects examples:",
			'  findObjects(area > 200, df) # identify objects with area > 200',
			'  findObjects(positive == TRUE, df, mask) # mark positive objects in mask')
		cat(usage, sep = "\n")
		return(invisible(NULL))
	}
# check 'df' argument
	if (!is(df, "data.frame"))
		stop("'", deparse(substitute(df)), "' must be a data.frame")
	dnames <- names(df)
	if (!("frame" %in% dnames & ("well" %in% dnames | "file" %in% dnames)))
		stop('"frame" and either "well" or "file" must be in ', deparse(substitute(df)))

# evaluate expression and assign to original data.frame
	sel <- eval(substitute(expr), df)
	vname <- tail(make.unique(c(names(df), "var")), 1)
	df[[vname]] <- sel

# split 'var' according to 'well/file' and 'frame' in 'df'
	if ("well" %in% names(df))
   group <- "well"
	else if ("file" %in% names(df))
		group <- "file"
	spl.1 <- split(df, df[[group]], drop = FALSE)
	spl.2 <- lapply(spl.1, function(v) split(v[[vname]], v[["frame"]], drop = TRUE))

	neg <- lapply(spl.2, function(v) lapply(v, function(x) which(x == FALSE)))
	pos <- lapply(spl.2, function(v) lapply(v, function(x) which(x == TRUE)))

# return answer if mask is absent
	if (is.null(mask))
		return(if(invert) neg else pos)
	
# helper function to check on compatibility between mask and data.frame 
	.Nobs <- function(mask)
	{
		dm <- dim(mask)
		if (length(dm) == 2)
			return(max(mask))
		else
			return(as.vector(apply(mask, 3, max)))		
	}

# check for same number of objects in mask and data.frame
	if (is(mask, "Image"))
		Nmask <- .Nobs(mask)
	else if(is(mask, "list"))
		Nmask <- as.vector(unlist(lapply(mask, .Nobs)))
	else
		stop("unable to use this combination of 'df' and 'mask'")

	Ndata <- as.vector(unlist(lapply(spl.2, lengths)))
	if (!identical(Nmask, Ndata))
		stop("'mask' and 'data' have different number of objects")

# process
	if (is(mask, "Image")) {
		pos <- unlist(pos, recursive = FALSE)
		neg <- unlist(neg, recursive = FALSE)
		ans <- if(invert) rmObjects(mask, pos) else rmObjects(mask, neg)
	}
	else if (is(mask, "list"))
		ans <- if(invert) Map(rmObjects, mask, pos) else Map(rmObjects, mask, neg)
	return(ans)
}
