#########################################################################################
# 
# Calculate the likely cutoff value by Otsu's method for a robust mean and sd for a 
# normal density in a mixture. If 'left', assumes background population lies to the left
# of the normal mean.
#
# returns single number 
#
#########################################################################################

find.bgnd <- function(x, mult=5, log=TRUE, left=TRUE)
{
	library(genefilter)
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
