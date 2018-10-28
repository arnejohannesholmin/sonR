#*********************************************
#*********************************************
#' Returns a list of the following strings: (1) the path to the event, (2) the event name, (3) the event number, (4) the path to the cruise, and (5) the cruise name.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param dir.type  is the name of the directory holding the data files (usually one of "tsd" and "raw")
#' @param ...  is used in agrep() for locating events based on approximate string matching.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD pathparts rm.na
#' @importFrom utils tail
#'
#' @export
#' @rdname event.path
#'
event.path <- function(event=1, cruise=2009116, esnm="MS70", dir.type="tsd", ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############### LOG: ###############
	# Start: 2012-07-07 - Clean version.

	
	##################################################
	##################################################
	########## Preparation ##########
	# Functions used for expanding the file path to the folders specifying the 'esnm' and the 'dir.type':
	getPresent_esnm = function(event, esnm){
		present_esnm = basename(list.dirs(event))
		if(length(present_esnm)==0){
			stop(paste("No directories in ", event,sep=""))
		}
		if(!any(tolower(present_esnm)==tolower(esnm))){
			warning(paste("The specified 'esnm' (",esnm,") is not present in ",event,". The first chosen: ",present_esnm[1],sep=""))
			file.path(event, present_esnm[1])
		}
		else{
			file.path(event, esnm)
		}
	}
	getPresent_dir.type = function(event, dir.type){
		present_dir.type = basename(list.dirs(event))
		if(length(present_dir.type)==0){
			stop(paste("No directories in ", event,sep=""))
		}
		if(!any(tolower(present_dir.type)==tolower(dir.type))){
			warning(paste("The specified 'dir.type' is not present in ",event,". The first chosen: ",present_dir.type[1],sep=""))
			file.path(event, present_dir.type[1])
		}
		else{
			file.path(event,dir.type)
		}
	}
	
	# 'event' may be given as a list:
	if(is.list(event)){
		event = event$event
	}
	event = event[1]
	
	# Somtimes it is simpler to only check that the event exists:
	#if(length(event)==1 && isTRUE(file.info(event)$isdir)){
	#	return(event)
	#}
	# The event may be given in a form not recognized by the function, in which case it is returned as it is, with a warning. For this the input is stored:
	NULLevent = FALSE
	if(length(event)==0){
		NULLevent = TRUE
		event = 1
	}
	if(suppressWarnings(is.na(as.numeric(event)))){
		warningevent = list(event=event)
	}
	else{
		warningevent = list()
	}
	

	########## Execution ##########
	# If the event ends with a folder called "tsd" or "raw", things are assumed to be ok:
	if(tolower(basename(as.character(event[1]))) %in% c("raw","tsd")){
		event = event[1]
	}
	
	if(is.na(file.info(event)$isdir)){
		if(!is.na(file.info(dirname(event))$isdir)){
			warning(paste("The folder \"",dir.type,"\" not found in the directory \"",dirname(event),"\""))
		}
		else{
			warning("Event \"",event,"\" not found")
		}
	}
	
	# Locate the folder "Events" and expand the path to reach the folder 'dir.type':
	event_subdir = pathparts(event)
	atEvent = suppressWarnings(max(which(tolower(substr(event_subdir, 1, 6))=="events")))
	
	# Expand to the expected folder:
	if(length(event_subdir)==atEvent){
		#warning("The event not suffuciently specified. The first event chosen")
		eventlist = list.files(event)
		event = file.path(event, eventlist[1])
		# Get the list of esnm:
		event = getPresent_esnm(event,esnm)
		# Get the list of data file types:
		event = getPresent_dir.type(event,dir.type)
	}
	else if(length(event_subdir)==atEvent+1){
		# Get the list of esnm:
		event = getPresent_esnm(event,esnm)
		# Get the list of data file types:
		event = getPresent_dir.type(event,dir.type)
	}
	else if(length(event_subdir)==atEvent+2){
		# Get the list of data file types:
		event = getPresent_dir.type(event,dir.type)
	}
	
	if(is.na(file.info(event)$isdir)){
		warning("Event \"",event,"\" not found. Maybe the 'dir.type' in misspelled")
	}
	
	# Strip the event of double "/":
	event = gsub("//","/",event)
	
	
	########## Output ##########
	# Split the event path into separate directory names:
	event_subdir = pathparts(event)
	eventname = event_subdir[atEvent+1]
	eventnr = NA
	cruisename = NA
	cruise = NA
	understroke = gregexpr("_",eventname)[[1]]
	atE = gregexpr("_E",eventname)[[1]][1]
	# Remove additional descriptional strings past the event number in the file name:
	if(any(is.na(understroke) | is.na(atE))){
		#warning("Events not listed in 'event'. Path not sufficiently accurate")
	}
	else{
		if(isTRUE(tail(understroke,1)>atE)){
			eventname = substr(eventname,1,min(understroke[understroke>atE])-1)
		}
		eventnr = strsplit(eventname,"")[[1]]
		suppressWarnings(isnotnumeric<-which(is.na(as.numeric(eventnr))))
		eventnr = paste(eventnr[seq(atE+2, min(isnotnumeric[isnotnumeric>atE+2]-1,nchar(eventname)))],collapse="")
		# Get the cruise path and cruise name:
		cruise = paste(event_subdir[seq_len(atEvent-1)],collapse=.Platform$file.sep)
		cruisename = event_subdir[max(which(tolower(event_subdir)=="events"),which(tolower(event_subdir)=="cases"))-1]
		understroke = gregexpr("_",cruisename)[[1]]
		if(length(understroke)==1){
			cruisename = substr(cruisename,1,understroke[1]-1)
		}
	}
	
	# If the 'event' was NULL, it was set to 1 to locate the cruise, and then it is set to NULL again here:
	out = list(event=event, eventname=eventname, eventnr=eventnr, cruise=cruise, cruisename=cruisename, esnm=basename(dirname(event)))
	
	if(length(cruise)==0){
		out["cruise"] = list(NULL)
		out["cruisename"] = list(NULL)
	}
	if(NULLevent){
		out["event"] = list(NULL)
		out["eventname"] = list(NULL)
		out["eventnr"] = list(NULL)
	}
	out
	##################################################
	##################################################
}
