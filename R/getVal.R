#' Get Value from ComputeFeatures.xxx Function
#' 
#' Extract a value from one of the \code{computeFeatures} family of
#' functions using a mask or a list of masks and target image
#' or list of images.  
#'
#' @param mask An object mask or list of objects masks with
#'   connected pixels having the same integer value. 
#' @param ref A fluorescent \code{Image} object or list of
#'   \code{Image}s corresponding to the objects in \code{mask}
#'   or \code{NULL} for \code{computeFeatures.shape}.
#' @param val A character string identifying the parameter
#'   to return from the computeFeatures function \code{FUN}. The
#'   default value of \code{"b.mean"} returns the mean intensity.
#' @param FUN A \code{\link[EBImage]{computeFeatures}} function to be
#'   applied over \code{mask} and \code{ref}. The default of 
#'   \code{NULL} uses the character string in \code{val} to identify
#'   the appropriate function.
#' @param simplify If \code{TRUE} (default), the result will be
#'   collapsed into a vector. Otherwise, a vector for each member of 
#'   the list \code{mask} and \code{ref} will be returned.
#' @param ... additional parameters are passed to \code{FUN}.
#' 
#' @details
#'
#' Objects identified in \code{mask} will be projected onto \code{ref} 
#' and quantified with the function specified in \code{FUN}. The single 
#' value specified by \code{val} will be returned as a list if \code{
#' simplify == FALSE} or as a single vector if \code{simplify == TRUE}.
#'
#' Common usages include getting the mean intensity from a \code{mask}
#' and reference image in \code{ref} by \code{getVal(mask, ref)} or
#' getting the area of objects from an object mask in \code{mask} by
#' \code{getVal(mask, val = "area")}.
#' 
#' Functions that can be explicitly specified include 
#' \code{\link[EBImage:computeFeatures]{computeFeatures.basic}}, 
#' \code{\link[EBImage:computeFeatures]{computeFeatures.shape}}, 
#' \code{\link[EBImage:computeFeatures]{computeFeatures.moment}} and 
#' \code{\link[EBImage:computeFeatures]{computeFeatures.haralick}}.
#' See the \code{EBImage} help page for appropriate values for
#' \code{FUN} and \code{val}. Note that "\code{basic quantiles}" 
#' parameters are misnamed with an extra '0'. For example, the 5th
#' percentile and 95th percentile values are named "\code{b.q005}"
#' and "\code{b.q095}".
#'
#' @return
#' 
#' A vector if \code{simplify} is \code{TRUE} or a list of values obtained
#' from \code{FUN}.
#' 
#' @import EBImage  
#' 
#' @export
#'
getVal <- function(mask, ref = NULL, val = "b.mean", FUN = NULL, simplify = TRUE, ...)
{
# current variables from computeFeatures functions
	varList <- list(
		basic = c("b.mad", "b.mean", "b.q001", "b.q005", "b.q05", "b.q095",
			"b.q099", "b.sd"),
		shape = c("s.area", "s.perimeter", "s.radius.max", "s.radius.mean",
			"s.radius.min", "s.radius.sd"),
		moment = c("m.cx", "m.cy", "m.eccentricity", "m.majoraxis", "m.theta"))

# help function if mask is missing
	if (missing(mask)) {
		cat("Accepted values for 'val' include:\n\n")
		print(varList)
		return(invisible(NULL))
	}

# assign computeFeatures function and variable
	if (is.null(FUN)) { # identify function from 'val'
		sel <- lapply(varList, function(v) grep(val, v, ignore.case = TRUE, value = TRUE))
		if (sum(lengths(sel)) == 0)
			stop("unable to match 'val' with computeFeatures function")
		if (sum(lengths(sel)) > 1)
			stop("'val' does not identify a unique computeFeatures variable")
		val <- unlist(sel)
		FUN <- switch(names(val),
			basic = EBImage::computeFeatures.basic,
			shape = EBImage::computeFeatures.shape,
			moment = EBImage::computeFeatures.moment)
	}
# assign 'ref' for case of computeFeatures.shape
	if (identical(FUN, EBImage::computeFeatures.shape)) {
		if (is(mask, "list"))
			ref <- rep(list(NULL), length(mask))
		else
			ref <- NULL
	}
# working function
	.proc <- function(mask, ref, val, FUN, ...) {
		if (is.null(ref))
			unlist(lapply(getFrames(mask), function(m) FUN(m, ...)[,val]))
		else
			unlist(Map(function(m, t, ...) FUN(m, t, ...)[,val], 
				getFrames(mask), getFrames(ref), ...))
	}
# determine if 'mask' and 'ref' are lists or Images
	if (is(mask, "Image") && (is(ref, "NULL") | is(ref, "Image")))
		ans <- .proc(mask, ref, val, FUN, ...)
	else if (is(mask, "list")) {
		if (!all(sapply(mask, is, "Image")))
			stop("'mask' must be an Image object or list of Image objects")
		if (length(mask) != length(ref))
			stop("'mask' and 'ref' must of the same length")
		ans <- Map(.proc, mask, ref, val, list(FUN)) # list() needed 
	}
	else 
			stop("unable to handle combination of 'mask' and 'ref'")
# unlist the results if simplify is TRUE
	if (simplify == TRUE)
		ans <- unlist(ans)
	return(ans)
}
