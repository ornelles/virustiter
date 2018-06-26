#########################################################################################
# 
# Calculate the likely cutoff value using Otsu's method in the package genefilter 
# for a normal density mix. left==TRUE assumes background population lies to the left
# of the normal mean. This function can be replaced with more sophisticated versions.
#
# returns single number 
#
#########################################################################################

findBgnd <- function(x, mult=5, log=TRUE, left=TRUE)
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
