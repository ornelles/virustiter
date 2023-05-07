#' Plot Titer Fit from GLM Model
#'
#' Prepare a diagnostic plot of viral titer GLM fit.
#'
#' @param fm Fitted model or list of fitted models from \code{\link{getFit}}.
#' @param line.col Color for best-fit line appropriate for \code{par("col")}.
#' @param main List or character vectors appropriate for the plot title
#'   of the same length as \code{fm}. If \code{NULL}, the current date will
#'   be used. 
#'
#' @details
#'
#' Base graphics are used to prepare the plot where the first argument is
#' either a single fitted model or a list of fitted models from
#' \code{\link{getFit}}. The function calls \code{\link{getEC63}} to obtain the fit
#' and confidence intervals.
#'
#' @return
#'
#' No value is returned. This function is called for the side-effect of 
#' producing a plot (or plots).
#'
#' @export
#'  
plotFit <-function(fm, line.col = 2, main = NULL)
{
	# ensure main is appropriate and pad with NULL if need be
		if (is.null(main))
			main <- list(NULL)
		main <- as.list(main)
		main <- c(main, rep(list(NULL), max(0, length(fm) - length(main))))
		main <- main[seq_along(fm)]
		
	# working function to dispatch on each glm model
	.plotFit <- function(fm, main)
	{
		if (is.null(main))
			main <- paste0("<", Sys.Date(), ">")

		moi <- exp(fm$model[[2]]) # model data.frame holds values used for fit
		y <- prop.table(fm$model[[1]],1)[,1]
		cf <- getEC63(fm)

		res <- fm$data	# entire set provided to glm()
		unit <- levels(res$unit)[1]
		txt <- sprintf("%0.3g %s (95%% CI:%0.3g-%0.3g) ", cf[1], unit, cf[2], cf[3])

		xlo <- min(moi[moi > 0])
		xhi <- max(moi)
		xp <- exp(seq(log(xlo), log(xhi), length = 101))
		yp <- predict(fm, data.frame(x = xp), type = "response")
		xpp <- cf[1]
		ypp <- 1 - exp(-1)

		plot(y ~ moi, subset = moi > 0, log = "x", las = 1, ylim = c(0, 1),
				xlab = paste("\n", "One IU = ", txt, sep = ""),
				ylab = "Infected fraction", main = main)
		lines(xp, yp, col = 2)
		lines(c(xlo, xpp, xpp), c(ypp, ypp, -0.02), lty = 2, col = 4)
	}
	if ("list" %in% class(fm))
		junk <- mapply(.plotFit, fm, main)
	else
		junk <- .plotFit(fm, main[[1]])
}
