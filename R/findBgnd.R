#' Find Background Value
#' 
#' Calculate a value between background and foreground.
#' 
#' @param x fluorescent values to evaluate
#' @param mult numeric multiplier applied to the standard deviation,
#'   with default value of 3
#' @param log a \code{logical} flag to use log-transformed values
#' @param crit the critical p-value between 0 and 1 to distinquish
#'   bimodal from non-bimodal populations, default of 0.1
#' 
#' @details
#' 
#' The \code{link[multimode]{modetest}} function provides a test of whether the
#' values are bimodal or multimodal. If the p-value from this test is less
#' than the critical value (\code{crit}), the background value is considered
#' to be the value at the first "valley" in the kernel density distribution.
#' For such bimodal distributions, the \code{mult} parameter is ignored.
#' 
#' If the distribution is not bimodal,the population is assumed to be a mixture
#' of values from a normal distribution of low values and a lognormal (or log)
#' distribution of high values. The maximum of low values is determined from a
#' kernel density estimate and the \emph{left} half of the distribution is fit
#' to a Gaussian distribution with the \code{\link[MASS]{fitdistr}} function.
#' The value returned is the position of the peak + \code{mult} times the
#' estimated standard deviation of the fitted distribution.
#' 
#' Because fluorescent values are typically log-transformed before analysis,
#' values are log-transformed by default. This can be turned
#' off with the \code{log} parameter. Typically, the parameter \code{mult}
#' must be empirically determined for non-bimodal distributions. (See the
#' examples for the impact of \code{mult} on the return value. Note that with
#' \code{mult = 0}, the mean value of the background intensity will be returned.
#' 
#' If the distribution is very heavily skewed to the left (mostly dark values),
#' the standard deviation of the distribution will be estimated as 1.35x
#' the interquartile range.
#'
#' It may be helpful to discard exceptionally low values before finding
#' an estimated background should the estimated background be absurdly low. 
#' 
#' @seealso \code{\link{getZero}}
#' @seealso \code{\link{getBgnd}}
#'
#' @return
#' 
#' For non-bimodally distributed values, the value estimated as position of the
#' most abundant values + mult * standard. For bimodally distributed values,
#' the value at the valley between the two populations.
#' 
#'  @examples
#'    dev.new(width = 6, height = 8)
#'    set.seed(1234)
#'    par(mfrow = c(2, 1))
#'    x1 <- c(rlnorm(100), rlnorm(30, 4))
#'    plot(density(log(x1)), main = "Bimodal ('mult' is ignored)")
#'    bg <- findBgnd(x1)
#'    abline(v = log(bg), col = 2) # background limit
#'    x2 <- c(rlnorm(80), rlnorm(40, 6, sdlog = 5))
#'    plot(density(log(x2)), main = "Non-bimodal")
#'    bg <- sapply(0:4, function(m) findBgnd(x2, mult = m))
#'    abline(v = log(bg), col = 1:5) # background limit
#'    txt <- expression(mult %*% sd)
#'    legend("topright", legend = 0:4, title = txt, lty = 1, col = 1:5)
#' 
#' @import
#' EBImage
#'
#' @importFrom MASS fitdistr
#' @importFrom multimode modetest locmodes
#'
#' @export
#' 
findBgnd <- function(x, mult = 3, log = TRUE, crit = 0.1, ratio.limit = 1/10)
{
	if (log) {
		if (all(x <= 0)) stop("positive values are needed if log = TRUE")
		x <- log(x[x > 0])
	}
	if (crit < 0 | crit > 1) stop("'crit' must be between 0 and 1")

# test for more than 1 peak (mode)
	bimodal <- FALSE
	pval <- multimode::modetest(x, mod0 = 1, method = "SI", B = 100)$p.value
	if (pval < crit) { # if less than crit, probably bimodal
		v <- suppressWarnings(multimode::locmodes(x, mod0 = 2))
		if (v$fvalue[1]/v$fvalue[3] > ratio.limit)
			bimodal <- TRUE
		else
			bimodal <- FALSE
	}

# find breakpoint
	if (bimodal == TRUE) 
		ans <- v$locations[2]
	else { # otherwise, see if left half can be fit to a Gaussian distribution
		d <- density(x[x > quantile(x, 0.01) & x < quantile(x, 0.99)]) # trim?
		xmid <- d$x[which.max(d$y)]
		xl <- x[x <= xmid]
		xx <- c(xl, 2*xmid  - xl) # extract left half
		if (length(xx) >= 0.05 * floor(length(x))) {
			fit <- MASS::fitdistr(xx, "normal")
			ans <- fit$est[1] + mult * fit$est[2]
		}
		else # if not, make best guess
			ans <- min(x) + mult * IQR(x, na.rm = TRUE)/1.349
	}

# return estimate
	if (log == TRUE)
		ans <- exp(ans)
	return(unname(ans))
}
