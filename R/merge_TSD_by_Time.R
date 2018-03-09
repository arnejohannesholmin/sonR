#*********************************************
#*********************************************
#' Merges the last file of the current file number and the first file of the next.
#'
#' @param files	list of file names, seed the use in \code{\link{compr.event}}.
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
merge_TSD_by_Time <- function(files, fun=weighted.mean, msg=FALSE, ...){
	
	# See info.TSD("nmtc").
	mergeUsingNmtc = function(name, fun, data, nmtc){
		temp = data[[name]]
		if(name %in% timevar){
			apply(temp, seq_len(length(dim(temp))-1), head, n=1)
			}
		else{
			apply(temp, seq_len(length(dim(temp))-1), fun, w=if(length(nmtc)==0) ones(tail(dim_all(temp)),1) else nmtc, ...)
			}
		}
		
	files = files[file.exists(files)]
	timevar = labl.TSD("t")
	# Read the first time step of each file:
	utim <- read.TSDs(files, var="utim", t=1, clean=FALSE)
	utim = unlist(utim[names(utim)=="utim"], use.names=FALSE)
	# Get the files with duplicated time steps:
	duputim = unique(utim[duplicated(utim)])
	common = which(utim %in% duputim)
	outfiles = files[!duplicated(utim)]
	commonfiles = files[common]
	commonutim = utim[common]
	
	if(msg){
		infostring = paste("Merging joint time steps:")
		cat(infostring,"\n",sep="")
		totalsteps = length(duputim)
		stepfact = nchar(infostring)/totalsteps
		oldvalue = 0
		}
	
	for(i in seq_along(duputim)){
		# Print a dot if the floor of the new value exceeds the old value:
		if(msg){
			thisvalue = floor(i*stepfact)
			if(thisvalue > oldvalue){
				cat(rep(".",thisvalue-oldvalue),if(i==totalsteps) "\n", sep="")
				oldvalue = thisvalue
				}
			}
		
		ind = which(commonutim == duputim[i])
		thisdata = read.TSDs(commonfiles[ind], drop.out=FALSE, clean=FALSE, header=FALSE)
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
		write.TSD(thisdata, commonfiles[ind][1], numt=1, ...)
		unlink(commonfiles[ind][-1])
		}
	outfiles
	}