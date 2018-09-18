#' Plot Titer Fit from GLM Model
#'
#' Prepare a diagnositc plot by column, row, or none.
#'
#' @param fm Fitted model or list of fitted models from \code{getFit()}.
#' @param by Character vector indicating that fits are organized as a
#'   single fit ("none") or by "column" or "row".
#' @param line.col Color for best-fit line (appropriate for \code{par("col")}.
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @details
#'
#' Base graphics are used to prepare the plot where the first argument is
#' either a single fitted model or a list of fitted models from
#' \code{getFit()}. The function calls \code{getTiter()} to obtain the fit
#' and confidence intervals.
#'
#' @return
#'
#' No value is returned. This function is called for the side-effect of 
#' producing a plot (or plots).
#'
#' @export
#'  
plotFit <-function(fm, by = c("none", "column", "row"), line.col = 2, ...)
{
	.plotFit <- function(fm) {				# internal function to work on glm fitted model
		moi <- exp(fm$model[[2]])				# model data.frame holds values used for fit
		y <- prop.table(fm$model[[1]],1)[,1]
		cf <- getTiter(fm)
		info <- try(sapply(well.info(rownames(fm$model)), unique), silent = TRUE)
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

		if (!("directory" %in% names(fm$data)))
			main <- paste(Sys.Date(), main.text)
		else
			main <- paste(fm$data$directory[1], main.text)

		res <- fm$data	# entire data.frame handed to glm()
		unit <- levels(res$unit)[1]
		txt <- sprintf("%0.3g %s (95%% CI:%0.3g-%0.3g) ", cf[1], unit, cf[2], cf[3])

		xlo <- with(res, min(moi[moi > 0]))
		xhi <- with(res, max(moi))
		xp <- exp(seq(log(xlo), log(xhi), length = 101))
		yp <- predict(fm, data.frame(moi = xp), type = "response")
		xpp <- cf[1]
		ypp <- 1 - exp(-1)
		plot(y ~ moi, subset = moi > 0, log = "x", las = 1, ylim = c(0,1),
				xlab = paste("\n", "One IU = ", txt, sep = ""),
				ylab = "Infected fraction", main = main, ...)
		lines(xp, yp, col = 2)
		lines(c(xlo, xpp, xpp),c(ypp, ypp, -0.02), lty = 2, col = 4)
	}
	if ("list" %in% class(fm))
		junk <- lapply(fm, .plotFit)
	else
		junk <- .plotFit(fm)
}
