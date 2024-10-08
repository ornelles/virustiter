% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFit.R
\name{plotFit}
\alias{plotFit}
\title{Plot Titer Fit from GLM Model}
\usage{
plotFit(
  fm,
  main = NULL,
  pch = 1,
  col = 1,
  col.pch = col,
  lty.fit = 1,
  col.fit = col,
  lty.ref = 2,
  col.ref = "gray",
  xlim = NULL,
  ylim = NULL,
  ann = par("ann"),
  axes = TRUE,
  xlab = NULL,
  ylab = NULL,
  frame.plot = axes,
  digits = 2,
  ...
)
}
\arguments{
\item{fm}{fitted model or list of fitted models from \code{\link{getFit}}.}

\item{main}{plot title as a character, vector of characters or list of
characters to be used as the plot title(s). Values will be replicated as
necessary. If \code{NULL}, the plot date will be used.}

\item{pch, col, col.pch}{plot character and color for data points. The default
color for the points and the fitted line is \code{col} but these can be
changed with the appropriate argument. \code{col} will also be used as the
default color for other plot elements.}

\item{lty.fit, col.fit}{line type and color for the GLM best-fit line.  Set
\code{lty.fit} to \code{NA} to exclude the best-fit line.}

\item{lty.ref, col.ref}{line type and color for the value on the x-axis
that intersect the 63% value on the best-fit line. Set \code{lty.ref} to
\code{NA} to exclude the reference line.}

\item{xlim, ylim}{optional parameters to override the default parameters of
the full range of x values for \code{xlim} and \code{c(0, 1)} for
\code{ylim}.}

\item{ann}{a logical value indicating whether the default annotation (title
and x and y axis labels) should appear on the plot.}

\item{axes}{a logical value indicating whether both axes should be drawn on
the plot.}

\item{xlab}{a label for the x axis. If the \code{unit} attribute associated
with \code{fm} is a volume (ml, ul, nl, pl or fl), the default will be the
titer as IU per ml. Otherwise, the default will be the x-value
corresponding to an MOI of 1 (63\% point).}

\item{ylab}{a label for the y axis, defaults to "Infected fraction".}

\item{frame.plot}{a logical value indicating whether a box should be drawn
around the plot.}

\item{digits}{integer indicating the number of significant decimal places 
to use for the calculated x-axis label}

\item{...}{arguments to be passed to \code{\link[graphics]{plot.default}}}
}
\value{
No value is returned. This function is called for the side-effect of 
producing a plot (or plots).
}
\description{
Prepare an informative plot of viral titer from the GLM fit.
}
\details{
Base graphics are used to prepare a plot with \code{log = "x"} from a single 
fitted model or a list of fitted models generated by \code{\link{getFit}}. 
The function calls \code{\link{getTiter}} to derive the titer and includes 
this information on the plot by default. The graphic parameters can be 
used to change or exclude many elements of the plot. See the help 
information in \code{\link[graphics]{plot.default}} for details.

If \code{fm} is a list of fitted models, each will be plotted. Use 
\code{par(ask = TRUE)} to see each in turn or use 
\code{par(mfrow = c(nr, nc))} to place \code{nr x nc} plots on one device.

Additional curves can be added with the \code{\link{addOneFit}} as shown in 
the examples.
}
\examples{
# sample result from tally() with two sets of data
  plotFit_data <- read.csv(system.file("extdata", "plotFit_data.csv", 
    package = "virustiter"))
  fm <- getFit(plotFit_data, by = "row")

# Default format and example with minimal annotation
  plotFit(fm$A)
  plotFit(fm$A, ann = FALSE, axes = FALSE, lty.ref = NA)

# Example of two plots
  plotFit(fm$A, "Two fits", xlab = "Multiplicity (ml/cell)")
  addOneFit(fm$B)
  legend("topleft", legend = c("A", "B"), lty = 1, col = c(1, 2), bty = "n")

# With two panel figure 
  par(mfrow = c(2, 1))
  plotFit(fm, main = names(fm))

}
