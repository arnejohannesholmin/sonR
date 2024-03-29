% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanSv.TSD.R
\name{meanSv.TSD}
\alias{meanSv.TSD}
\title{Locates the point at which the kernel desity estimate of the Sv of the subset is equal to half its maximum, and uses Paper III of the PhD of Holmin to estimate the mean Sv from this. The output variables all start with X, denoting the SX90 segmentation method, but are applicable to all segmentation methods.}
\usage{
meanSv.TSD(
  data,
  plot.hist = FALSE,
  allow.vbsc = TRUE,
  list.out = FALSE,
  minlen = 1,
  type = "H",
  enlarged = FALSE,
  ...
)
}
\arguments{
\item{data}{A list containing the acoustic data eihter given as 'vbss' for a segment of the data or as 'vbsc'.}

\item{plot.hist}{Logical: If TRUE plot the histogram of the Sv.}

\item{allow.vbsc}{?}

\item{list.out}{Use a list as output.}

\item{minlen}{The minimum length of the data, at and below which NA is returned.}

\item{type}{A single character denoting the type of segmentation to label the output with.}

\item{enlarged}{Logical: If TRUE add "E" to the variable names and remove the last character.}

\item{...}{Passed on to \code{\link{density}} and \code{\link{hist}}.}
}
\description{
Locates the point at which the kernel desity estimate of the Sv of the subset is equal to half its maximum, and uses Paper III of the PhD of Holmin to estimate the mean Sv from this. The output variables all start with X, denoting the SX90 segmentation method, but are applicable to all segmentation methods.
}
\examples{
\dontrun{}

}
