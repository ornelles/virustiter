#' Find Cutoff by Otsu's Method
#' 
#' Calculate a working cutoff value using Otsu's method in the package
#' \code{genefilter} for a mixture of values with a normal density. The
#' method assumes the background population lies to the left
#' of the normal mean of the population. This can be changed by the parameter
#' \code{left}.
#' 
#' @param x Fluorescent values to evaluate.
#' @param mult Arbitrary numeric multiplier applied to standard deviation.
#' @param log \code{logical} flag to use log-transformed values.
#' @param left \code{logical} flag to indicate that the background population lies to the
#' left (\code{TRUE}) of the distribution.
#' 
#' @details
#' 
#' Data are assumed to be drawn from two unimodal, continuous 
#' distributions where the background population of a normal distribution 
#' lies on the left. The mode of this population is estimated by an 
#' iterative, half-range method. The cutoff value returned is the mode +
#' \code{mult} * standard deviation. Note that fluorescent values are
#' typically log-transformed before analysis. The \code{mult} must be
#' empirically determined. The default value of 5 is fairly conservative. 
#' 
#' @return
#' 
#' The cutoff estimate as half.mode + mult * standard deviation of the normal population.
#' 
#' @import
#' EBImage
#' genefilter
#' 
#' @examples
#'   x <- c(rnorm(1000), rnorm(500, mean = 3))
#'   findBgnd(x, mult = 0, log = FALSE)
#'   plot(density(x))
#'   abline(v = findBgnd(x, mult = 2, log = FALSE))
#'
#' @export
#' 
findBgnd <- function(x, mult = 5, log = TRUE, left = TRUE)
{
	if(require(genefilter) == FALSE)
		stop("load 'genefilter' package or replace 'nucMask' function")
	if (log) {
		x <- x[x>0]
		x <- log(x)
	}
	x.mu <- half.range.mode(x)
	if (left) {
		side <- (x - x.mu)[x < x.mu]
		x.sd <- sqrt(sum(side^2)/(length(side)-1))
		ret <- x.mu + mult*x.sd
	}
	else {
		side <- (x - x.mu)[x > x.mu]
		x.sd <- sqrt(sum(side^2)/(length(side)-1))
		ret <- x.mu - mult*x.sd
	}
	if (log)
		ret <- exp(ret)
	return(ret)
}
