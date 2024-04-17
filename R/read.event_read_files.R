#*********************************************
#*********************************************
#' Reads all files of a specific type. Used in read.event().
#'
#' @param files  is a vector of the file names of the files.
#' @param filesind  is a vector of the indexes of the files in the list of files.
#' @param tlist  is the list of time indexes to be read for each file.
#' @param var  is a vector of the variables to read.
#' @param origin  is (1) a vector of two elements representing the origin of the global coordinate system (G), (2) the numbering index of the ping in the total sequence of pings of the event, which is to be regarded as the origin of (G) (ignoring heave so that the x-y-plane of (G) is on the surface of the sea), or (3) NULL, implying that the origin be put to the mid point of the vessel posistions. Used only for vessel files.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD global2car read.TSDs
#'
#' @export
#' @rdname read.event_read_vessel
#'
read.event_read_files<-function(files, filesind, tlist, var, origin=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-08-01 - Clean version.
	# Update: 2013-08-09 - Added origin==NULL, implying that the origin be put to the mid point of the vessel posistions.
	# Last: 2016-07-07 - Created from the old read.event_read_vessel().
	########### DESCRIPTION: ###########
	# Reads all vessel-files. Used in read.event().
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	

	##################################################
	##################################################
	# Read the vessel information:
	# merge=TRUE combines vessel information from multiple files into vectors, and indt=FALSE is a consequence of the way 'tlist' is defined (read.event() always considers general time step indexes and 'tlist' is defined according to this requirement):
	out = suppressWarnings(read.TSDs(files[filesind], t=tlist[filesind], var=var, merge=TRUE, indt=FALSE, msg=FALSE, drop.out=FALSE))
	
	# Only for vessel files:
	# Calculate the cartesian positions of the vessel in (G) defined by 'origin':
    #browser()
	psxvRequested = ("all" %in% var && file_ext(files[filesind][1])=="vessel") || any(c("psxv", "psyv") %in% var)
	psxvNotPresent = any(c(length(out$psxv)==0, length(out$psyv)==0))
	if(psxvRequested && psxvNotPresent){
		temp = suppressWarnings(read.TSDs(files[filesind], t=tlist[filesind], var=c("lonv","latv"), merge=TRUE, indt=FALSE, msg=FALSE, drop.out=FALSE))[c("lonv","latv")]
		out[names(temp)] = temp
		if(length(origin)==1){
			# clean=TRUE preserves only the first occurence of the longitude and latitude information, and indt=TRUE is used because read.event() considers general time step indexes:
			origin = suppressWarnings(read.TSDs(files[filesind], t=origin, var=c("lonv","latv"), clean=TRUE, indt=FALSE, msg=FALSE))
			}
		#if(length(origin)>0){
			out$psyv = global2car(cbind(c(out$lonv), c(out$latv)), origin)
			out$psxv = out$psyv[, 1]
			out$psyv = out$psyv[, 2]
		#	}
		}
	
	out
	##################################################
	##################################################
	}
