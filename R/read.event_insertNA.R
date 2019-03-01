#*********************************************
#*********************************************
#' Adds empty list elements of NAs (if vectors) to data read in read.event(). This is done by comparing to the time step indices stored in the object 'TIME'.
#'
#' @param thisout  is a list containing the relevant data.
#' @param TIME  is the list returned from UNIX_time().
#' @param filesind  is a vector of the indexes of the relevant files in the list of files.
#' @param filelist  is the list of files.
#' @param t  is the tine step indexes to return data at.
#' @param tlist  is the list of time step indexes for the files.
#' @param var  is a vector of the variables to to insert NAs to.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD NAs rm.na
#'
#' @export
#' @rdname read.event_insertNA
#'
read.event_insertNA <- function(thisout, TIME, filesind, filelist, t, tlist, var){
	
	############### LOG: ###############
	# Start: 2013-05-12 - Clean version.
	# Last: 2013-10-01 - Removed 'merge'.
	# Get the specific time step indexes of each file:
	INDT <- vector("list", length(TIME$indt[filesind]))
	for(i in seq_along(INDT)){
		INDT[[i]] <- TIME$indt[filesind][[i]][tlist[filesind][[i]]]
		}
	INDT <- unique(unlist(INDT))
	
	# The information about which file the variables returned by read.TSDs() are located in is given in the optional field 'nvarFile' (remove this to avoid trouble with the following code):
	nvarFile <- thisout$nvarFile
	thisout$nvarFile <- NULL
	# Get the number of time steps in each of the variables, returned to 'thisnumt':
	thisnumt <- sapply(thisout, function(x) if(length(dim(x)>1)) ncol(x) else length(x))
	# Check whether the length of any of the variables in 'thisout' are shorter than the number of requested time steps:
	insertNA <- which(thisnumt>1 & thisnumt<length(t) & (names(thisout) %in% var))
	
	# Move through the variables which are too short, and insert the data at the correct time steps:
	if(length(insertNA)>0){
		for(i in insertNA){
			# Get the number of the file from which the variable was read:
			if(is.list(thisout[[i]])){
				temp <- vector("list", length(t))
				temp[rm.na(match(INDT,t))] <- thisout[[i]]
				thisout[[i]] <- temp
				}
			else if(length(dim(thisout[[i]]))==0){
				temp <- NAs(length(t))
				temp[rm.na(match(INDT,t))] <- thisout[[i]]
				thisout[[i]] <- temp
				}
			}
		}
	thisout
}
