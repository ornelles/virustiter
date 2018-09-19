#
# Uses half population at left to estimate normal (log-normal) distribution
# and then infer the mean and sd. Rather like half-range mode of 'genefilter'
#
# I think I like this one better and like the idea of using the multiplier
# as a multiplier for the estimated SD. More logical than "5"...
#
findBgnd2 <- function(x, k = 3, log = TRUE)
{
	if (log == TRUE) x <- log(x[x > 0])
	d <- density(x)
	xmid <- d$x[which.max(d$y)]
	xl <- x[x <= xmid]
	xx <- c(xl, 2*xmid  - xl)
	fit <- MASS::fitdistr(xx, "normal")
	ans <- fit$est[1] + k * fit$est[2]
	if (log == TRUE) ans <- exp(ans)
	return(ans)
}