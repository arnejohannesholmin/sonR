#*********************************************
#*********************************************
#' Sorts in a case insensitive matter. Created to solve the problem of errors when listing files in R as run from the terminal, in which case sorting is not case insensitive. In the R-GUI sorting is done case insensitive.
#'
#' @param x  is an R-object as input to sort() / order().
#' @param decreasing  used in sort() / order().
#' @param na.last  used in sort() / order().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname list.files_caseInsensitive
#'
list.files_caseInsensitive<-function(path, pattern=NULL, all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE, method="shell"){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-10-12 - Clean version.
	########### DESCRIPTION: ###########
	# Sorts in a case insensitive matter. Created to solve the problem of errors when listing files in R as run from the terminal, in which case sorting is not case insensitive. In the R-GUI sorting is done case insensitive.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is an R-object as input to sort() / order().
	# ---decreasing--- used in sort() / order().
	# ---na.last--- used in sort() / order().
	
	
	##################################################
	##################################################
	path = list.files(path=path, pattern=pattern, all.files=all.files, full.names=full.names, recursive=recursive, ignore.case=ignore.case)
	# Set the locale to "en_US.UTF-8" before listing files:
	old_CTYPE = Sys.getlocale("LC_CTYPE")
	old_COLLATE = Sys.getlocale("LC_COLLATE")
	
	#newEnc <- "en_US.UTF-8"
	newEnc <- "C"
	if(old_CTYPE!=newEnc){
		Sys.setlocale("LC_CTYPE", newEnc)
		}
	if(old_COLLATE!=newEnc){
		Sys.setlocale("LC_COLLATE", newEnc)
		}
	# Order:
	path = path[order(tolower(path), method=method)]
	# Restore the previous locale:
	Sys.setlocale("LC_CTYPE", old_CTYPE)
	Sys.setlocale("LC_COLLATE", old_COLLATE)
	path
	##################################################
	##################################################
	}
