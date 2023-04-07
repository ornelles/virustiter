#' Get cutoff values for nuclear mask area
#'
#' Find appropriate range limits for the area of nuclei defined in a
#' binary nuclear mask. The limits will be defined by the "saddle"
#' point on the low end and by the upper quantile on the high end.
#'
#' @param mask An \code{Image} object, array, or a \code{list} of
#'   these objects containing a binary mask defining regions of the image
#'   to be segmented. Typically, this would obtained from
#'   \code{\link[virustiter]{nucMask}}. 
#' @param lo A optional numeric value indicating the approximate minimum value.
#'   This value will be used to search for the true saddle point. If negative,
#'   the search will start at the lower 10 percentile of the area.
#' @param hi A optional numeric value for the expected upper limit of the
#'   cutoff. If missing, the\code{quantile} will be used to determine the
#'   upper limit.
#' @param breaks Number of bins used for the area histograms (default of 256).
#' @param quant If \code{hi} is not provided, the area given by this
#'   quantile will be used as the upper limit.
#' @param plot If \code{TRUE}, plot histogram with selected cutoff points
#'
#' @details
#'
#' Get the low and high cutoffs to give to \code{\link[virustiter]{trimMask}}.
#'
#' @return
#' 
#' Vector of two values recommended as nuclear mask cutoff parameters.
#'
#' @import EBImage
#'
#' @export
#'
	getCutoff <- function(mask, lo, hi, breaks = 256, quant = 0.995, plot = FALSE)
	{
		area <- getVal(mask, val = "area")
		if (missing(hi)) hi <- quantile(area, quant)
		if (missing(lo) || lo < 0) lo <- quantile(area, 0.1) # use 10% as a guess
		hh <- hist(area[area < quantile(area, quant)], breaks = breaks, plot = FALSE)
		adj <- min(which(hh$mids > lo))
		j <- which.max(hh$counts[hh$mids > lo]) + adj
		i <- which.min(runmed(hh$counts[1:j], 5))
		ans <- c(hh$mids[i], hi)
		if (plot) {
			hist(area, breaks = breaks)
			abline(v = ans, col = 2)
		}
		return(unname(ans))
	}
