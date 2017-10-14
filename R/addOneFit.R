#########################################################################################
#
# add single fitted model with base graphics-provides options to add labels
#
#########################################################################################

addOneFit <- function(fm, line.col=4, ref.col=4, pch.col=line.col, ...) {
	moi <- exp(fm$model[[2]])		# model data.frame holds values used for fit
	y <- prop.table(fm$model[[1]],1)[,1]
	res <- fm$data  # entire data.frame handed to glm()
	cf <- getTiter(fm)

	xlo <- with(res, min(moi[moi > 0]))
	xhi <- with(res, max(moi))
	xp <- exp(seq(log(xlo), log(xhi), length=101))
	yp <- predict(fm, data.frame(moi=xp), type="response")
	xpp <- cf[1]
	ypp <- 1-exp(-1)

	points(y ~ moi, subset=moi>0, col=pch.col, ...)
	lines(xp, yp, col=line.col)
	lines(c(xlo,xpp,xpp),c(ypp,ypp,-0.02), lty=2, col=ref.col)
}
