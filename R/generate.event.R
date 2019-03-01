#*********************************************
#*********************************************
#' Generate a simulation event.
#'
#' @param event  is the name of the event, such as "S2014119_D200141023_E0003".
#' @param cruise  is the name of the cruise, such as "S2014119_PG.O.Sars[4174]".
#' @param esnm  is a vector of four character upper case names of the acoustic instruments, such as "MS70".
#' @param dir.type  is the types of data, usually c("raw", "tsd").
#' @param ...  Passed on to dir.create.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname generate.event
#'
generate.event<-function(event="Event1", cruise="Cruise1", esnm=NULL, dir.type=c("raw", "tsd"), ...){
	
	############### LOG: ###############
	# Start: 2014-02-21 - Clean version.
	
	##### Preparation #####
	if(!file.exists(cruise)){
		cruise <- Acoustics_datasets_directory()
	}
	#if(is.null(dir.data)){
	#	dir.data=Acoustics_datasets_directory()
	#	}
	
	
	##### Execution and output #####
	if(length(cruise>0)){
		#event = c(outer(file.path(dir.data, cruise, "Events", event, esnm), dir.type, file.path))
		event = c(outer(file.path(cruise, "Events", event, esnm), dir.type, file.path))
		}
	suppressWarnings(lapply(event, dir.create, recursive=TRUE, ...))
	event
}
