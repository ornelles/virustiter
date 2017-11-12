#
# p2p interactively measure point to point distance on an image (plot)
#
p2p <- function(n = 512, col = "magenta", type = "o", pch = 3, ...)
{
	ans <- numeric()
	while (n > 0) {
		p <- locator(2, type = type, pch = pch, col = col, ...)
		if (is.null(p))
			break
		ans <- c(ans, sqrt(sum(sapply(p, diff)^2)))
		n <- n - 1
	}
	return(ans)
}
