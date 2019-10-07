#' Get Value from ComputeFeatures.xxx Function
#' 
#' Extract a value through one of the \code{computeFeatures} family of
#' functions using an object mask and reference image or a list of each.
#'
#' @param mask An object mask or list of objects masks with
#'   connected pixels having the same integer value. 
#' @param ref A fluorescent \code{Image} object or list of
#'   \code{Image} objects corresponding to the objects in \code{mask}
#'   or \code{NULL} if a reference image is not required.
#' @param val A character string identifying the parameter
#'   to return from the \code{computeFeatures} function assigned
#'   to \code{FUN}. The default value of \code{"b.mean"} returns the
#'   mean intensity with the
#'   \code{\link[EBImage:computeFeatures]{computeFeatures.basic}}
#'   function.
#' @param FUN A \code{\link[EBImage]{computeFeatures}} function to be
#'   applied over \code{mask} and \code{ref}. The default of 
#'   \code{NULL} uses the character string in \code{val} to select
#'   the appropriate function. If this function is specified,
#'   the value in the argument \code{val} must be an exact match
#'   to one of the values returned by the specified function.
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
#' Common usages include getting the mean intensity from an object mask
#' and reference image  or getting the area of objects from an object mask
#' as shown below.
#' \preformatted{
#' # To extract the mean intensity:
#'   mfi <- getVal(mask, ref) # uses default parameter "b.mean"
#'
#' # To extract object area, collapsed to a single vector:
#'   area <- getVal(mask, val = "area") # no need 'ref'
#'
#' # A complex example of Haralick difference entropy at
#' # at a Haralick scale of 4 pixels
#'   den <- getVal(mask, ref, val = "h.den.s4",
#'     FUN = computeFeatures.haralick, haralick.scales = 4)
#' }
#' 
#' The character string in \code{val} is used as a \code{grep} pattern to
#' search among the parameters returned by the \code{computeFeatures}
#' functions. Consequently, a partial match such as "ecc" is sufficient
#' for a parameter such as "m.eccentricity". However, because "mean" occurs
#' in "b.mean" and "s.radius.mean", a partial match to "mean" does not
#' identify a unique parameter. The unique parameter is used to assign the
#' value to \code{FUN}. This value can be explicitly assigned as one of
#' \code{\link[EBImage:computeFeatures]{computeFeatures.basic}}, 
#' \code{\link[EBImage:computeFeatures]{computeFeatures.shape}}, 
#' \code{\link[EBImage:computeFeatures]{computeFeatures.moment}} and 
#' \code{\link[EBImage:computeFeatures]{computeFeatures.haralick}}.
#' See the \code{EBImage} help page for appropriate values for
#' \code{FUN} and \code{val}. Note that "\code{basic quantiles}" 
#' parameters are mislabeled with an extra '0'. For example, the 5th
#' percentile and 95th percentile values are named "\code{b.q005}"
#' and "\code{b.q095}", respectively.
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
# capture additional arguments in ...
	dots <- list(...)
	if (length(dots) == 0) dots <- list(NULL)

# error check on 'val'
	if(!(is(val, "character") && length(val) == 1))
		stop("'val' must be a single character string")

# assign computeFeatures function and variable
	if (is.null(FUN)) { # identify function from 'val'
		sel <- lapply(varList, function(v) grep(val, v, ignore.case = TRUE, value = TRUE))
		if (sum(lengths(sel)) == 0)
			stop("unable to match ", val, " with a computeFeatures function")
		if (sum(lengths(sel)) > 1)
			stop(val, " does not identify a unique computeFeatures variable")
		val <- unlist(sel)
		FUN <- switch(names(val),
			basic = EBImage::computeFeatures.basic,
			shape = EBImage::computeFeatures.shape,
			moment = EBImage::computeFeatures.moment)
	}
# assign 'ref' to 'NULL' for case of computeFeatures.shape
	if (identical(FUN, EBImage::computeFeatures.shape)) {
		if (is(mask, "list"))
			ref <- rep(list(NULL), length(mask))
		else
			ref <- NULL
	}
# working function
	.proc <- function(mask, ref, val, FUN, dots) {
		if (is.null(ref)) {
			res <- lapply(getFrames(mask), function(m)
				do.call(FUN, args = c(list(m), dots))[,val])
			return(unlist(res))
		}
		else {
			res <- Map(FUN, getFrames(mask), getFrames(ref), MoreArgs = dots)
			res <- do.call(rbind, res)
			return(res[,val])
		}
	}
# determine if 'mask' and 'ref' are lists or Images
	if (is(mask, "Image") && (is(ref, "NULL") | is(ref, "Image")))
		ans <- .proc(mask, ref, val, FUN, dots)
	else if (is(mask, "list")) {
		if (!all(sapply(mask, is, "Image")))
			stop("'mask' must be an Image object or list of Image objects")
		if (length(mask) != length(ref))
			stop("'mask' and 'ref' must of the same length")
		ans <- Map(.proc, mask, ref, val, list(FUN), list(dots)) # list() needed 
	}
	else 
			stop("unable to handle combination of 'mask' and 'ref'")
# unlist the results if simplify is TRUE
	if (simplify == TRUE)
		ans <- unlist(ans)
	return(ans)
}
