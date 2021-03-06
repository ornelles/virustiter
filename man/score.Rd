% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{score}
\alias{score}
\title{Score Cells as Positive/Negative}
\usage{
score(df, bgnd = NULL, param = "mfi", ...)
}
\arguments{
\item{df}{A data.frame generated by \code{\link{parseImages}} or
and an appropriate grouping variable.}

\item{bgnd}{A numeric vector of length one or a named numeric vector for
each level of the grouping variable, typically \code{well, row, column, or file}.
If missing, the function \code{\link{getBgnd}} will be called with parameters
provided in \code{...} to determine the background.}

\item{param}{A character string identifying the variable to be scored.}

\item{...}{Additional arguments for \code{\link{getBgnd}} such as \code{by}.}
}
\description{
Assign \code{logical} variable \code{positive} as \code{TRUE} or
\code{FALSE} in a data.frame of fluorescent intensity values based
on the cutoff values in the argument \code{bgnd}.
}
\details{
A variable named \code{positive} will be added to the data.frame if 
necessary and assigned \code{TRUE} or \code{FALSE} if the value 
named in \code{"param"} is greater than the appropriate value in
\code{bgnd}. If \code{bgnd == NULL}, \code{\link{getBgnd}} will be called
with any additional arguments in \code{...} to determine the background.
}
\examples{
# Small set of data by stack
  f <- system.file("extdata", "by_stack", package = "virustiter")
  i <- getImages(f, pattern = "file005")
  v <- parseImages(i)

# Score with grouping value of "file"
  v <- score(v, by = "file")
  xtabs(~ positive, v)
  plot(log(mfi) ~ area, v, col = c(1,2)[positive + 1])

}
