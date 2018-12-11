#' Find Background Value by Otsu's Method
#' 
#' Calculate a working cutoff value using an approximation of Otsu's method
#' for a mixture of values with a normal density. This function assumes the
#' background population has a normal distribution and lies to the left of
#' distribution. 
#' 
#' @param x Fluorescent values to evaluate.
#' @param mult Numeric multiplier applied to standard deviation.
#' @param log \code{logical} flag to use log-transformed values.
#' 
#' @details
#' 
#' Data are assumed to be drawn from two unimodal, continuous 
#' distributions where the background population follows a normal distribution. 
#' The maximum of this population is determined from a kernel density estimate
#' and the \emph{left} half of the distribution is fit to a Gaussian
#' distribution. The value returned is the position of the peak +
#' \code{mult} times the standard deviation of the distribution.
#' Note that fluorescent values are typically log-transformed before analysis.
#' The parameter \code{mult} must be empirically determined.
#'
#' This code replaces a previous version that used the \code{half.range.mode}
#' function in the \code{genefilter} package.
#' 
#' @return
#' 
#' The value estimated as position of the background peak + mult * standard
#' deviation of the normal population.
#' 
#' @import
#' EBImage
#'
#' @importFrom MASS fitdistr
#' 
#' @examples
#'   x <- c(rlnorm(200), rlnorm(120, 4))
#'   plot(density(log(x)))
#'   abline(v = log(findBgnd(x, mult = 0)), col = 4) # peak
#'   abline(v = log(findBgnd(x)), col = 2) # background limit
#'
#' @export
#' 
findBgnd <- function(x, mult = 2.5, log = TRUE)
{
	if (log == TRUE)
		x <- log(x[x > 0])
	d <- density(x)
	xmid <- d$x[which.max(d$y)]
	xl <- x[x <= xmid]
	xx <- c(xl, 2*xmid  - xl)
	fit <- MASS::fitdistr(xx, "normal")
	ans <- fit$est[1] + mult * fit$est[2]
	if (log == TRUE)
		ans <- exp(ans)
	return(ans)
}
