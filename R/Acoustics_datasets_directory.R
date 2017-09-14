#*********************************************
#*********************************************
#' Set the directory of the acoustic data used by the sonR package
#'
#' @param path	The path to the directory in which to put acoustic data, structured as CruiseName/"Events"/EventName/EchosounderName/"raw" for the raw files and substituting "raw" by "tsd" for TSD files.
#'
#' @export
#'
Acoustics_datasets_directory <- function(path=NULL){
	file = "Acoustics_datasets_directory.txt"
	extdata = file.path(find.package("sonR"), "extdata")
	file = file.path(extdata, file)
	if(length(path)>0 && nchar(path)>0){
		return(writeLines(path, file))
	}
	extdatafiles = list.files(extdata, full.names=TRUE)
	
	out <- NULL
	warnmessage <- "Use \n\nAcoustics_datasets_directory(PATH_TO_DIRECTORY_HOLDING_ACOUSTIC_DATA_READ_BY_THE_PACKAGE_sonR) \n\nto set the path to the directory holding resource files and acoustic data from cruises (such as \"~/Data/Acoustics\")"
	if(file %in% extdatafiles){
		out <- readLines(file, warn=FALSE)
	}
	else{
		stop(paste0("The file \n\n", file, "\n\nmissing. ", warnmessage))
	}
	if(length(out)==0){
		warning(paste0("The file ", file, " is empty. ", warnmessage))
	}
	return(out)
}
