#*********************************************
#*********************************************
#' Sorts in a case insensitive matter. Created to solve the problem of errors when listing files in R as run from the terminal, in which case sorting is not case insensitive. In the R-GUI sorting is done case insensitive.
#'
#' @param ...		Used in \code{\link{list.files}}.
#' @param method	The ordering method used in \code{\link{order}}.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname list.files_caseInsensitive
#'
list.files_caseInsensitive <- function(path=".", pattern=NULL, all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE, include.dirs=FALSE, no..=FALSE, method="shell"){
	
	############### LOG: ###############
	# Start: 2012-10-12 - Clean version.
	# Last: 2018-11-17 - Cleaned up.
	
	path <- list.files(path=path, pattern=pattern, all.files=all.files, full.names=full.names, recursive=recursive, ignore.case=ignore.case, include.dirs=include.dirs, no..=no..)
	# Set the locale to "en_US.UTF-8" before listing files:
	old_CTYPE <- Sys.getlocale("LC_CTYPE")
	old_COLLATE <- Sys.getlocale("LC_COLLATE")
	
	#newEnc <- "en_US.UTF-8"
	newEnc <- "C"
	if(old_CTYPE != newEnc){
		Sys.setlocale("LC_CTYPE", newEnc)
		}
	if(old_COLLATE != newEnc){
		Sys.setlocale("LC_COLLATE", newEnc)
		}
	# Order:
	path <- path[order(tolower(path), method=method)]
	# Restore the previous locale:
	Sys.setlocale("LC_CTYPE", old_CTYPE)
	Sys.setlocale("LC_COLLATE", old_COLLATE)
	path
}
