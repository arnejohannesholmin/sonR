#*********************************************
#*********************************************
#' This function returns a list of indices of files of common time (utim)/sailed distance (sadv), or if this list is given merges the common files.
#'
#' @param files	list of file names, seed the use in \code{\link{compr.event}}.
#' @param var  The variable used as the merging variable, e.g., "utim" or "sadv".
#' @param ind  A list of indices of common files. Each element of the list contains two or more file indices of files to be merged.
#' @param fun  funciton specifying how to merge the information in the common time interval.	
#' @param msg  Logical: if TRUE print status bar.
#' @param ...  Passed on to \code{fun} and \code{\link{write.TSD}}.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom utils tail head
#' @importFrom stats weighted.mean
#'
#' @export
#'
merge_TSD_by_var <- function(files, var="utim", ind=NULL, fun=weighted.mean, msg=FALSE, ...){
	
	# See info.TSD("nmtc").
	mergeUsingNmtc <- function(name, fun, data, nmtc){
		
		timevar <- labl.TSD("t")
		
		temp = data[[name]]
		if(name %in% timevar){
			apply(temp, seq_len(length(dim(temp))-1), head, n=1)
		}
		else{
			apply(temp, seq_len(length(dim(temp))-1), fun, w=if(length(nmtc)==0) ones(tail(dim_all(temp)),1) else nmtc, ...)
		}
	}
	mergeOneGroup <- function(ind, files, fun){
		thisdata = read.TSDs(files[ind], drop.out=FALSE, clean=FALSE, header=FALSE)
		# Check that all variables have the same number of list elements (identically structured files):
		tableOfNames = table(names(thisdata))
		unames = setdiff(unique(names(thisdata)), "nmtc")
		if(!all(tableOfNames==tableOfNames[1])){
			warning("All files must have the same variables")
			}
		nmtc = unlist(thisdata[names(thisdata)=="nmtc"], use.names=FALSE)
		thisdata = lapply(unames, function(x) mergeListKeepDimensions(thisdata[names(thisdata)==x]))
		names(thisdata) = unames
		thisdata = lapply(unames, mergeUsingNmtc, fun=fun, data=thisdata, nmtc=nmtc)
		names(thisdata) = unames
		thisdata$nmtc = sum(nmtc)
		write.TSD(thisdata, files[ind][1], numt=1, ...)
		unlink(files[ind][-1])
	}
	
	# Assure that the files exist:
	files <- files[file.exists(files)]
	
	# Get and return the indices of common files:
	if(length(ind)==0){
		# Read the first time step of each file:
		utim <- read.TSDs(files, var=var, t=1, clean=FALSE)
		utim <- unlist(utim[names(utim)==var], use.names=FALSE)
		# Get the files with duplicated time steps:
		duputim <- unique(utim[duplicated(utim)])
		common <- which(utim %in% duputim)
		commonfiles <- files[common]
		commonutim <- utim[common]
	
		# Get the file indices of common times:
		ind <- lapply(duputim, function(this) which(utim == this))
		ind
	}
	# If ind is already given, merge the files:
	else{
		temp <- papply(ind, mergeOneGroup, files=files, fun=fun, msg=msg, pb=FALSE)
		toRemove <- unlist(lapply(ind, "[", -1))
		outfiles <- files[-toRemove]
		outfiles
	}
}
