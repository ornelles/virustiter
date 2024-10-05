#' Find Background Value
#' 
#' Find the fluorescent value between negative and positive cells
#' 
#' @param x fluorescent values to evaluate
#' @param mult numeric multiplier applied to the standard deviation
#'   for a non-bimodal distributions, default value of 3
#' @param log a \code{logical} flag to use log-transformed values
#' @param crit the critical p-value between 0 and 1 to distinguish
#'   bimodal from non-bimodal populations, default of 0.1
#' @param ratio.limit if the ratio of the first peak height to the second peak
#'   height exceeds this value (1/10), the population could be bimodal
#' 
#' @details
#' 
#' The \code{\link[multimode]{modetest}} function tests whether the values have
#' a bimodal (or multimodal) distribution. If the p-value from this test is less
#' than the critical value (\code{crit}), the background value is considered
#' to be the value at the first "valley" in the kernel density distribution.
#' For such bimodal distributions, the \code{mult} parameter is ignored.
#' 
#' If the distribution is not bimodal, the population is assumed to be a mixture
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
#' must be empirically determined for non-bimodal distributions. See the
#' examples for the impact of \code{mult} on the return value. Note that with
#' \code{mult = 0}, the peak value of a non-bimodal distribution will be
#' returned.
#' 
#' If the distribution is very heavily skewed to the left (mostly dark values),
#' the standard deviation of the distribution will be estimated as 1.35x
#' the interquartile range.
#'
#' It may be helpful to discard exceptionally low values before finding
#' an estimated background should the estimated background be absurdly low.
#' The code currently discards the upper and lower 1 per cent of values before
#' fitting to a normal (or lognormal) distribution.  
#' 
#' @seealso \code{\link{getZero}}
#' @seealso \code{\link{getBgnd}}
#'
#' @return
#' 
#' For bimodally distributed values, the value at the valley between the two
#' populations. For non-bimodally distributed values, the value of the
#' most abundant value + mult * standard deviation of the estimated distribution.
#' 
#' @import
#' EBImage
#'
#' @importFrom MASS fitdistr
#' @importFrom multimode modetest locmodes
#'
#' @examples
#'   dev.new(width = 6, height = 8)
#'   set.seed(1234)
#'   par(mfrow = c(2, 1))
#'
#'   x1 <- c(rlnorm(100), rlnorm(30, 4))
#'   bg <- findBgnd(x1)
#'   plot(density(log(x1)), main = "Bimodal ('mult' is ignored)")
#'   abline(v = log(bg), col = 2) # background limit
#'
#'   x2 <- c(rlnorm(80), rlnorm(40, 6, sdlog = 5))
#'   bg <- sapply(0:4, function(m) findBgnd(x2, mult = m))
#'   plot(density(log(x2)), main = "Non-bimodal")
#'   abline(v = log(bg), col = 1:5) # background limit
#'   txt <- expression(mult %*% sd)
#'   legend("topright", legend = 0:4, title = txt, lty = 1, col = 1:5)
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

# Is there evidence for a biomodal distribution?
# This only tests for the likelihood of more than 1 peak (mode).
# Distributions with more than two peaks would be treated as non-bimodal.
# The method of Silverman (method = "SI") was chosen with 100 replicates
# for the sake of speed.
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

# return estimated breakpoint
  if (log == TRUE)
    ans <- exp(ans)
  return(unname(ans))
}
