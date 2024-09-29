#' Get cutoff values for nuclear mask area
#'
#' Find appropriate range limits for the area of nuclei defined in a
#' binary nuclear mask. The limits will be defined by the "saddle"
#' point on the low end and by the upper quantile on the high end.
#'
#' @param area either an \code{Image} object, array, or a \code{list} of
#'   these objects containing a binary mask defining regions of the image
#'   to be segmented. Typically, these binary masks would obtained from
#'   \code{\link[virustiter]{nucMask}}. Alternatively, this can be a numeric
#'   value for which "cutoff" values are desired.
#' @param lo a optional numeric value indicating the approximate minimum value.
#'   This value will be used to search for the true saddle point. If negative,
#'   the search will start at the lower 10 percentile of the area.
#' @param hi a optional numeric value for the expected upper limit of the
#'   cutoff. If missing, the\code{quantile} will be used to determine the
#'   upper limit.
#' @param breaks a number (or means) to determining the number of cells for
#'   the optional histogram (see \code{link[graphics]{hist}})
#' @param quant if \code{hi} is not provided, the area given by this
#'   quantile will be used as the upper limit.
#' @param plot if \code{TRUE}, plot histogram with selected cutoff points
#'
#' @details
#'
#' Get the low and high cutoffs to give to \code{\link[virustiter]{trimMask}}.
#'
#' @return
#' 
#' Vector of two values recommended as cutoff values suitable for
#' \code{link[virustiter]{trimMask}}
#'
#' @import EBImage
#'
#' @export
#'
	getCutoff <- function(area, lo, hi, breaks = 256, quant = 0.995, plot = FALSE)
	{
	# argument check
		if (is(area, "list")) {
			sel <- sapply(area, function(x) is.integer(imageData(x)))
			if (!all(sel))
				stop("'", deparse(substitute(area)),
					"' is a list but not all are integer Image masks")
			else
				area <- getVal(area, val = "s.area")
		}
		else if (is(area, "Image") && is.integer(imageData(area))) {
			if (!is.integer(imageData(area)))
				stop("'", deparse(substitute(x)), "' is not an integer Image mask")
			else
				area <- getVal(area, val = "s.area")
		}
		else if (!is.numeric(area))
			stop("unable to use '", deparse(substitute(x)), "'")

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
