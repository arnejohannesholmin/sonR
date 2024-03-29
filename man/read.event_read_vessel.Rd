% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.event_read_files.R
\name{read.event_read_files}
\alias{read.event_read_files}
\title{Reads all files of a specific type. Used in read.event().}
\usage{
read.event_read_files(files, filesind, tlist, var, origin = NULL)
}
\arguments{
\item{files}{is a vector of the file names of the files.}

\item{filesind}{is a vector of the indexes of the files in the list of files.}

\item{tlist}{is the list of time indexes to be read for each file.}

\item{var}{is a vector of the variables to read.}

\item{origin}{is (1) a vector of two elements representing the origin of the global coordinate system (G), (2) the numbering index of the ping in the total sequence of pings of the event, which is to be regarded as the origin of (G) (ignoring heave so that the x-y-plane of (G) is on the surface of the sea), or (3) NULL, implying that the origin be put to the mid point of the vessel posistions. Used only for vessel files.}
}
\description{
Reads all files of a specific type. Used in read.event().
}
\examples{
\dontrun{}

}
