% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.event.R
\name{generate.event}
\alias{generate.event}
\title{Generate a simulation event.}
\usage{
generate.event(event = "Event1", cruise = "Cruise1", esnm = NULL,
  dir.type = c("raw", "tsd"), dir.data = NULL, ...)
}
\arguments{
\item{event}{is the name of the event, such as "S2014119_D200141023_E0003".}

\item{cruise}{is the name of the cruise, such as "S2014119_PG.O.Sars[4174]".}

\item{esnm}{is a vector of four character upper case names of the acoustic instruments, such as "MS70".}

\item{dir.type}{is the types of data, usually c("raw", "tsd").}

\item{schoolthr}{is a function of size and mean volume backscattering strength (S_V) defining the threshold level below the 90-percentile of the school using the initial above-noise-threshold. Can also be a vector of values, where one segmentation file is written for each value of 'schoolthr'.}
}
\value{

}
\description{
Generate a simulation event.
}
\examples{
\dontrun{}

}