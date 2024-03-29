% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernSmooth1.R
\name{kernSmooth1}
\alias{kernSmooth1}
\title{(Usually Gaussian) smoothing an array along the first dimension. Note: Due to high variability in the speed of the fft() used in convolve() in this function, the performance times vary greatly for different lengths of the input data. NAs are removed in the actual convolution, so a large number of NAs should not slow down the function.}
\usage{
kernSmooth1(x, kern = 3, nsd = 3)
}
\arguments{
\item{x}{is an array to be smoothed along the first dimension.}

\item{kern}{is the kernel of the smoothing filter, defaulted to a Gaussian kernel.}

\item{nsd}{is the number of standard deviations on either side of the mean of the Gaussian kernel, outside which the kernel is zero.}
}
\description{
(Usually Gaussian) smoothing an array along the first dimension. Note: Due to high variability in the speed of the fft() used in convolve() in this function, the performance times vary greatly for different lengths of the input data. NAs are removed in the actual convolution, so a large number of NAs should not slow down the function.
}
\examples{
\dontrun{}

}
