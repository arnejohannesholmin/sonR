% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ksmooth.SG.TSD.R
\name{ksmooth.SG.TSD}
\alias{ksmooth.SG.TSD}
\title{Applies a Gaussian 1-D, 2-D, or 3-D kernel on the 'data' given as a list with names according to the TSD-convension, and the bandwidths given in the vector 'h':}
\usage{
ksmooth.SG.TSD(
  data,
  h = 1,
  w = h * 3,
  sim = FALSE,
  ind = list(),
  drop = FALSE,
  var = "vbsc",
  na.rm = 1
)
}
\arguments{
\item{data}{is the list of data with names according to the TSD-convension.}

\item{h}{is the three element vector of bandwidths.}

\item{w}{is the boundary of the kernel, outside which the it is 0.}

\item{sim}{is a TRUE if smoothing should be done only along the first dimensions, simultaneously over the stages of the last dimension. If 'sim' is an integer larger than 1, the positions 'coords' are used 'sim' times, and the data 'x' should have length 'sim' times the length of one coordinate of 'coords'.}

\item{ind}{is a list of indexes, as given to subset_TSD(), used to select the subset over which the estimation of high intensity noise is done. Defaulted to exclude the first 100 voxels along each beam.}

\item{drop}{is TRUE if the data should be dropped of unused dimensions prior to smoothing. This can be useful for smoothing data which .}

\item{var}{is a string specifying the variable to smooth.}

\item{na.rm}{is single integer representing the dimension along which NAs are discarded from the smoothing in the case that sim==1. For example, if na.rm=2 and the dimension of 'x' is [5,12,7], and x[3,2:4,5]=NA, then all data x[,2:4,] will be excluded from the smoothing and set to NA. If na.rm=FALSE, no NAs should be contained in the data.}
}
\description{
Applies a Gaussian 1-D, 2-D, or 3-D kernel on the 'data' given as a list with names according to the TSD-convension, and the bandwidths given in the vector 'h':
}
\examples{
\dontrun{}

}
