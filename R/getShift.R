#' Determine Optimal Image Shift
#' 
#' Determine the optimal x-y shift to align the target image with the
#' (nuclear) mask. 
#' 
#' @param mask An object mask with connected pixels having the same 
#'   integer value. 
#' @param target A fluorescent \code{Image} object corresponding to the 
#'   nuclear masks in \code{mask}. 
#' @param parscale A numeric vector of length 2 used as the 
#'   \code{parscale} parameter for the \code{\link{optim}} function. See 
#'   \code{\link{optim}} for more details. 
#' @param FUN Function to be minimized over \code{x[1]} and \code{x[2]} 
#'   by \code{\link[stats]{optim}}. See the details for more information.  
#' 
#' @details
#' 
#' This optimization function \code{FUN} accepts three arguments 
#' (\code{x, mask, target}) where \code{x} is a numeric vector of
#' length two representing the x and y position in the integer 
#' \code{Image} mask (\code{mask}) and the \code{Image} object 
#' (\code{target}). \code{FUN} must return a single value
#' representing the difference between the two images. The 
#' \code{\link[stats]{optim}} function minimizes \code{FUN} with respect to
#' \code{x[1]} and \code{x[2]}. The default function is designed
#' to align nuclear masks with predominantly nuclear signals and
#' should be replaced for other localization patterns. 
#' 
#' The argument \code{parscale} is used by \code{\link[stats]{optim}} to scale the 
#' parameters in \code{x} such that a unit change in the parameter 
#' amounts to unit change in the optimizing function empirically, 
#' \code{c(25, 25)} seems to be in the middle of a robust range.  
#' 
#' @return
#' 
#' A list of subpixel translations named \code{"dx"} and \code{"dy"} 
#' that can be applied to the argument \code{target} with
#' \code{\link{translate}} to maximize the alignment between mask and target. 
#' 
#' @examples
#'  path <- system.file("extdata", "by_folder/b2", package = "virustiter")
#'  x <- getImages(path)
#'  getShift(nucMask(x$nuc[[1]]), x$tgt[[1]])
#'
#' @import EBImage  
#' 
#' @export
#' 
getShift <- function(mask, target, parscale = c(25, 25), FUN = idiff)
{
	dm <- dim(mask)
	if (!identical(dm, dim(target)))
		stop("mask and target are of different sizes")
	if (length(parscale) == 1)
		parscale <- rep(parscale, 2)

# working function
	.getShift <- function(mask, target, FUN, parscale) {
		res <- optim(c(0, 0), fn = FUN, mask = mask, target = target,
			control = list(parscale = parscale))
		setNames(round(res$par, 1), c("dx", "dy"))
	}

# return a list for an 2 x n array of shifts
	if (length(dm) > 2)
		ans <- lapply(seq_len(dm[3]),
			function(i) .getShift(mask[,,i], target[,,i], FUN, parscale))
	else 
		ans <- as.list(.getShift(mask, target, FUN, parscale))
	return(ans)
}

#
# local function for optimization
#
idiff <- function(x, mask, target) {
	dm <- dim(mask)
	if (!identical(dm, dim(target)))
		stop("mask and target are of different sizes")
	if (length(dm) != 2)
		stop("idiff is meant to be called on single images")

	xm <- mask > 0 # convert integer mask to binary values
	xn <- normalize(target) # normalize target image between 0 and 1
	xt <- translate(xn, x, filter = "none") # translated normal image
	xp <- xt * !xm # include only those pixels outside of the mask
	if (all(xp == 0)) # just in case identical images were aligned
		return(0)
	else
		return(mean(xp[xp > 0])) # mean used rather than sum to avoid edge artifacts
}
