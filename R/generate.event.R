#*********************************************
#*********************************************
#' Generate a simulation event.
#'
#' @param event  is the name of the event, such as "S2014119_D200141023_E0003".
#' @param cruise  is the name of the cruise, such as "S2014119_PG.O.Sars[4174]".
#' @param esnm  is a vector of four character upper case names of the acoustic instruments, such as "MS70".
#' @param dir.type  is the types of data, usually c("raw", "tsd").
#' @param schoolthr  is a function of size and mean volume backscattering strength (S_V) defining the threshold level below the 90-percentile of the school using the initial above-noise-threshold. Can also be a vector of values, where one segmentation file is written for each value of 'schoolthr'. 
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname generate.event
#'
generate.event<-function(event="Event1", cruise="Cruise1", esnm=NULL, dir.type=c("raw", "tsd"), dir.data=NULL, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-02-21 - Clean version.
	########### DESCRIPTION: ###########
	# Segments SX90 data using noise-up and top-down.
	########## DEPENDENCIES: ###########
	# 
	############ VARIABLES: ############
	# ---event--- is the name of the event, such as "S2014119_D200141023_E0003".
	# ---cruise--- is the name of the cruise, such as "S2014119_PG.O.Sars[4174]".
	# ---esnm--- is a vector of four character upper case names of the acoustic instruments, such as "MS70".
	# ---dir.type--- is the types of data, usually c("raw", "tsd").
	# ---schoolthr--- is a function of size and mean volume backscattering strength (S_V) defining the threshold level below the 90-percentile of the school using the initial above-noise-threshold. Can also be a vector of values, where one segmentation file is written for each value of 'schoolthr'. 
	
	
	##################################################
	##################################################
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
	##################################################
	##################################################
	}
