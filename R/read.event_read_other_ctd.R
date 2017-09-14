#*********************************************
#*********************************************
#' Reads ctd-variables. Used in read.event().
#'
#' @param var  is a vector of the variables to read.
#' @param ctdfiles  is a vector of the file names of the ctd-files.
#' @param ctdfilesind is a vector of indices of the ctd-files in the list of files associated with 'TIME'.
#' @param TIME  is the list returned from UNIX_time().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSDs
#'
#' @export
#' @rdname read.event_read_other_ctd
#'
read.event_read_other_ctd<-function(var, ctdfiles, ctdfilesind, TIME){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-08-01 - Clean version.
	########### DESCRIPTION: ###########
	# Reads ctd-variables. Used in read.event().
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---var--- is a vector of the variables to read.
	# ---ctdfiles--- is a vector of the file names of the ctd-files.
	# ---TIME--- is the list returned from UNIX_time().
	

	##################################################
	##################################################
	suppressWarnings(read.TSDs(ctdfiles[sapply(TIME$labl[ctdfilesind],function(x) any(var %in% x))], var=var, clean=TRUE, msg=FALSE, t="all"))
	##################################################
	##################################################
	}
