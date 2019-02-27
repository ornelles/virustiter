#' Get Value from ComputeFeatures.xxx Function
#' 
#' Extract a value from one of the 'computeFeatures' family of
#' functions using a mask or a list of masks and target image
#' or list of images.  
#'
#' @param mask An object mask or list of objects masks with
#'   connected pixels having the same integer value. 
#' @param target A fluorescent \code{Image} object or list of
#'   \code{Image}s corresponding to the objects in \code{mask}. 
#' @param val A character string identifying the parameter
#'   to return from the computeFeatures function \code{FUN}. The
#'   default value of \code{"b.mean"} returns the mean intensity.
#' @param FUN A \code{\link[EBImage]{computeFeatures}} function to be
#'   applied over \code{mask} and \code{target}, default of 
#'   \code{computeFeatures.basic}.
#' @param simplify If \code{TRUE} (default), the result will be
#'   collapsed into a vector. Otherwise, a vector for each member of 
#'   the list \code{mask} and \code{target} will be returned. 
#' 
#' @details
#'
#' Objects identified in \code{mask} will be projected onto \code{target} 
#' and quantified with the function specified in \code{FUN}. The single 
#' value specified by \code{val} will be returned as a list if \code{
#' simplify == FALSE} or as a single vector if \code{simplify == TRUE}.
#' 
#' Possible functions include \code{computeFeatures.basic, 
#' computeFeatures.shape, computeFeatures.moment}, and \code{
#' computeFeatures.haralick}. See the help page for
#' \code{\link[EBImage]{computeFeatures}} for possible values for
#' \code{FUN} and \code{val}.
#'
#' @return
#' 
#' A vector if \code{simplify} is\code{TRUE} or a list of values obtained
#' from \code{FUN}.  
#' 
#' @import EBImage  
#' 
#' @export
#'
	getVal <- function(mask, target, val = "b.mean",
			FUN = computeFeatures.basic, simplify = TRUE)
	{
	# working function
		.proc <- function(mask, target, val, FUN) {
			dm <- dim(mask)
			if (length(dm) == 2) dim(mask) <- dim(target) <- c(dm, 1)
			unlist(lapply(seq_len(dim(mask)[3]),
					function(i) FUN(mask[,,i], target[,,i])[,val]))
		}
	# determine if 'mask' and 'target' are lists or Images
		if (is(mask, "Image") && is(target, "Image"))
			ans <- .proc(mask, target, val, FUN)
		else if (is(mask, "list") && is(target, "list")) {
			if (!all(sapply(mask, is, "Image")) || !all(sapply(target, is, "Image")))
				stop("'mask' and 'target' must be Image objects or lists of Image objects")
			if (length(mask) != length(target))
				stop("'mask' and 'target' must of the same length")
			if (!identical(sapply(mask, dim), sapply(target, dim)))
				stop("'mask' and 'target' have different dimensions")
			ans <- Map(.proc, mask, target, val, list(FUN)) # list() needed 
		}
		else 
				stop("unable to handle combination of 'mask' and 'target'")
	# unlist the results if simplify is TRUE
		if (simplify == TRUE)
			ans <- unlist(ans)
		return(ans)
	}

