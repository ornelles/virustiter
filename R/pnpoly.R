#
# point inclusion in polygon test (by W. Randolph Franklin)
#
# points - points to be tested as a form appropriate for xy.coords 
# vertices - points that define the vertices of the polygon,
#	also processed by xy.coords where order defines the polygon and the
#	first and last points are connected
#
pnpoly <- function(points, vertices)
{
	pp <- xy.coords(points)
	vv <- xy.coords(vertices)
	nvert <- length(vv$x)

  # working function
	.fun <- function(x, y, vx, vy, nvert) {
		inside <- FALSE
		j <- nvert
		for (i in seq_len(nvert)) {
			if (((vy[i] > y) != (vy[j] > y)) &&
					(x < (vx[j] - vx[i]) * (y - vy[i]) / (vy[j] - vy[i]) + vx[i]))
				inside <- !inside
			j <- i
		}
		return(inside)
	}

  # apply to each pair of coordinates in points
	sapply(seq_len(lengths(points)[1]),
		function(k) .fun(pp$x[k], pp$y[k], vv$x, vv$y, nvert))
}