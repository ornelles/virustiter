#
# adaptive threshold using disc instead of box
#
thresh2 <- function(x, width, offset, boundary = c("replicate", "circular")) {
	boundary <- match.arg(boundary)
	r <- width - width%%2 + 1
	f <- makeBrush(r, shape = "disc")
	f <- f/sum(f)
	return (x > (filter2(x, f, boundary) + offset))
}
