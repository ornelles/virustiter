#' Plot One Titer Fit with more Options
#'
#' Prepare a single diagnostic plot with more labeling options than \code{plotFit()}.
#'
#' @param fm Fitted model or list of fitted models from \code{getFit()}.
#' @param main,xlab,ylab,pch.col Parameters passed to \code{plot}.
#' @param line.col Color for best-fit line (appropriate for \code{par("col")}.
#' param ref.col Color for reference lines identifying key Poisson points.
#' @param by Character vector indicating that fits are organized as a
#'   single fit ("none") or by "column" or "row".
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @details
#'
#' Base graphics are used to prepare the plot where the first argument is
#' a single fitted model from \code{getFit()}. The function calls
#' \code{getTiter()} to obtain the fit and confidence intervals.
#'
#' @return
#'
#' No value is returned. This function is called for the side-effect of 
#' producing a plot with base graphics.
#'
#' @export
#'  
plotOneFit <- function(fm, main = NULL, xlab = NULL, ylab = NULL,
				pch.col = "black", line.col = "red", ref.col = "blue", ...)
{
	moi <- exp(fm$model[[2]])      # model data.frame holds values used for fit
	y <- prop.table(fm$model[[1]],1)[,1]
	cf <- getTiter(fm)
	info <- try(sapply(well.info(rownames(fm$model)),unique), silent = TRUE)
	if (class(info) != "try-error") {
		info.len <- sapply(info, length)
		if(info.len["row"] == 1)
			main.text <- paste("Row", info$row)
		else if (info.len["column"] == 1)
			main.text <- paste("Column", info$column)
		else
			main.text <- ""
	}
	else
		main.text <- ""
	if (is.null(main)){
		if (!("directory" %in% names(fm$data)))
				main <- paste(Sys.Date(), main.text)
		else
				main <- paste(fm$data$directory[1], main.text)
	}

	res <- fm$data  # entire data.frame handed to glm()
	unit <- levels(res$unit)[1]
	txt <- sprintf("%0.3g %s (95%% CI:%0.3g-%0.3g) ", cf[1], unit, cf[2], cf[3])

	xlo <- min(moi[moi > 0])
	xhi <- max(moi)
	xp <- exp(seq(log(xlo), log(xhi), length = 101))
	yp <- predict(fm, data.frame(x = xp), type = "response")
	xpp <- cf[1]
	ypp <- 1-exp(-1)
	if (is.null(xlab)) xlab <- paste("\n", "One IU = ", txt, sep = "")
	if (is.null(ylab)) ylab <- "Infected fraction"

	plot(y ~ moi, subset = moi > 0, log = "x", las = 1, ylim = c(0,1),
		xlab = xlab, ylab = ylab, main = main, col = pch.col, ...)
	lines(xp, yp, col = line.col)
	lines(c(xlo,xpp,xpp), c(ypp,ypp,-0.02), lty = 2, col = ref.col)
}
