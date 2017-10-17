#########################################################################################
# plotPlate
#
# display entire plate with lattice graphics (expects well values as "a03", etc.
#
#########################################################################################

plotPlate <- function(df, cex=1/2, alpha=1/2, main, invert.y=TRUE, ... ) {
	if (!"well" %in% names(df))
		stop("requires variable 'well' in ", deparse(substitute(df)))
	library(lattice)

	n <- nlevels(df$well)
	if (n > 96)							# 384-well plate
		{rows <- 32; columns <- 24}
	else if (n > 48)					# 96-well plate
		{rows <- 8; columns <- 12}
	else if (n > 24)					# 48-well plate
		{rows <- 6; columns <- 8}
	else if (n > 12)					# 24-well plate
		{rows <- 4; columns <- 6}
	else if (n > 6)						# 12-well plate
		{rows <- 3; columns <- 4}
	else if (n == 6)					# 6-well plate
		{rows <- 2; columns <- 3}
	else								# fewer than 6
		{rows <- 1; columns <- n}

	if (any(paste(letters[1:8],rep(2:4,each=8),sep="") %in% levels(df$well)))
		ww <- sprintf("%s%d", rep(letters[1:rows], each=columns), 1:columns)
	else
		ww <- sprintf("%s%02d", rep(letters[1:rows], each=columns), 1:columns)

	df$well <- factor(df$well, levels=ww)	# revise levels

	skip <- ! levels(df$well) %in% unique(as.character(df$well))

	if (missing(main) & !("directory" %in% names(df)))
		main <- Sys.Date()
	else
		main <- df$directory[1]

	if (!"positive" %in% names(df)) df$positive <- FALSE
	if (invert.y) df$ym <- -df$ym

	obj <- xyplot(ym ~ xm | well, data=df, group=positive, cex=cex, alpha=alpha,
			as.table=TRUE, aspect="iso", layout=c(columns,rows), skip=skip,
			xlab = "", ylab = "", scales = list(draw=FALSE), main = main, ...)
	plot(obj)
	invisible(obj)
}
