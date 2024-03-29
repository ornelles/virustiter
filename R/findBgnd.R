#' Find Background Value by Otsu's method or a best guess
#' 
#' Calculate a value between background and foreground pixels using
#' an approximation of Otsu's method for a bimodal mixture of values
#' or a unimodal normal distribution. This function assumes that the
#' background value will lie to the right of the left-most population.
#' 
#' @param x Fluorescent values to evaluate.
#' @param mult Numeric multiplier applied to the standard deviation,
#'   with default value of 2.5.
#' @param log A \code{logical} flag to use log-transformed values.
#' 
#' @details
#' 
#' Data are assumed to be drawn from two unimodal, continuous 
#' distributions where the background population follows a normal distribution. 
#' The maximum of this population is determined from a kernel density estimate
#' and the \emph{left} half of the distribution is fit to a Gaussian
#' distribution with the \code{\link[MASS]{fitdistr}} function. The value
#' returned is the position of the peak + \code{mult} times the estimated
#' standard deviation of the distribution.
#'
#' Because fluorescent values are typically log-transformed before analysis,
#' the values are log-transformed by default before analysis. This can be turned
#' off with the \code{log} parameter. Typically, the parameter \code{mult}
#' must be empirically determined. Note that with \code{mult = 0}, the mean
#' value of the background intensity will be returned.
#'
#' If the distribution is heavily skewed to the left (mostly dark values),
#' the standard deviation of the distribution will be estimated as 1.35x
#' the interquartile range. 
#'
#' @seealso \code{\link{getZero}}
#' @seealso \code{\link{getBgnd}}
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
#'   bg <- sapply(0:4, function(m) findBgnd(x, mult = m))
#'   abline(v = log(bg), col = 1:5) # background limit
#'   legend("topleft", legend = 0:4, title = "mult", lty = 1, col = 1:5)
#'
#' @export
#' 
findBgnd <- function(x, mult = 2.5, log = TRUE)
{
	if (log) {
		if (all(x <= 0)) stop("positive values are needed if log = TRUE")
		x <- log(x[x > 0])
	}
# determine density distribution and see if it is symmetrical
	d <- density(x)
	xmid <- d$x[which.max(d$y)]
	xl <- x[x <= xmid]
	xx <- c(xl, 2*xmid  - xl)

# fit to a half-distribution if at least 2% of the values are in the left
	if (length(xx) >= 0.02 * floor(length(x))) {
		fit <- MASS::fitdistr(xx, "normal")
		ans <- fit$est[1] + mult * fit$est[2]
	}

# otherwise use interquartile range assuming an exponential-like distribution
	else
		ans <- min(x) + mult * IQR(x, na.rm = TRUE)/1.349

# return estimate
	if (log == TRUE)
		ans <- exp(ans)
	return(unname(ans))
}
