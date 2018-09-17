#*********************************************
#*********************************************
#' Writes and returns the UNIX_time file of an event.
#'
#' @param event  is a string or vector of strings representing the path to the files or the directory from which to obtain the UNIX_time information, either by reading the UNIX_time file or by generating new information.
#' @param file  has one of three values: (1) The path to the file to which the UNIX time information should be written. (2) NULL/""/TRUE indicating the default file "UNIX_time.tsd". (3) FALSE, indicating that the file will not be written even though it would appear to be missing.
#' @param var  is a vector of the variable names of the variables to return. Available variables are "f000", "i000", "l000", "n000" and "r000".
#' @param t  is a vector of the time steps to return, which in this function corresponds to the files, as listed alphabetically (as returned from list.files(), with recursive=TURE)
#' @param fresh  is TRUE to force creating new UNIX time information.
#' @param allow.old  is a TRUE if old UNIX_time file is accepted (still with the correct list of files).
#' @param recursive  is FALSE to only consider files in the root tsd directory (used in files_caseInsensitive()).
#' @param saveData  is TRUE to save the UNIX_time info in memory for use later in the code (only recommended for events with a very high number of files).
#' @param cores  is the number of cores for parallel reading of the TSD files.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ones pathparts read.TSD strff write.TSD
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils tail head
#'
#' @export
#' @rdname UNIX_time
#'
UNIX_time <- function(event, file=FALSE, var="all", t="all", fresh=FALSE, msg=TRUE, allow.old=FALSE, recursive=TRUE, saveData=FALSE, cores=1){
	
	############### LOG: ###############
	# Start: 2011-08-11 - Clean version.
	# Update: 2012-07-31 - Changed to return full file names (not only base names, which inconveniently was the case before). Added the method of discarding beams-files and ctd-files when extracting the time point indexes 'i000', and setting i000=1 for these files. Also added the method of setting integer varibles to dtyp="shrt" (2 bytes) if the range is enclosed in +-32767. Finally added the file names to the variables 'u000', 'i000', 'l000' and 'n000'.
	# Update: 2012-08-03 - Fixed bug when returning utim, which was not discarded for beams-files and ctd-files.
	# Update: 2012-09-11 - Added the option 't'.
	# Update: 2012-10-12 - Changed to use the function list.files_caseInsensitive() instead of list.files(), due to the difference in the performance of the latter function in the terminal version of R and in the GUI version. The custom function list.files_caseInsensitive() sorts in a case insensitive matter, as is done by the GUI version of R.
	# Update: 2013-03-21 - Fixed bug where the UNIX_time.tsd file was unnecessarily re-written due to case dissimilarity between ald and new file names list. Also removed the info option.
	# Update: 2013-07-18 - Added long instead of floa to save integers.
	# Update: 2014-09-28 - Removed n000 and r000, and added t000.
	# Update: 2014-12-11 - Fixed some bugs and implemented that only the filea that have been changed are read.
	# Update: 2015-03-17 - Removed 'nottimefiles'.
	# Update: 2016-07-22 - Changed to compress sequences and remove u000. u000 is deduced from i000 and U000 when calling UNIX_time().
	# Last: 2016-11-03 - Fixed a major bug where new files were not registered in the event, and also set the default of 'saveData' to FALSE.
	
	
	##################################################
	##################################################
	##### Preparation #####
	processvar = c("f000", "u000", "i000", "l000", "t000", "U000", "I000", "i000_full", "I000_full")
	processdtyp = c("char", "doub", "long", "char", "char", "doub", "long", "long", "long")
	names(processdtyp) = processvar
	varNotWrite = c("u000", "i000_full", "I000_full")
	varToWrite = setdiff(processvar, varNotWrite)
	validvar = c("f000", "u000", "i000", "l000", "t000", "U000", "I000")
	validvarList = setdiff(validvar, c("U000", "I000"))
	if(identical(tolower(var),"all")){
		var = validvar
		}
	if(is.list(event)){
		event = event$event
		}
	if(length(event)==1 && isTRUE(file.info(event)$isdir)){
		dir = event
		}
	else{
		dir = dirname(event[1])
		}
		
	# Get the path to the UNIX_time file, if present:
	unixfile = file.path(dir,"UNIX_time.tsd")
	# If the file "UNIX_time.tsd" exists, extract the unix time points:
	if(!file.exists(unixfile)){
		unixfile = NULL
		}
	file_path_understroke = gsub("/", "_", file_path_sans_ext(unixfile))
	UNIX_time_data = paste0("UNIX_time_data", file_path_understroke)
	UNIX_time_mtime = paste0("UNIX_time_mtime", file_path_understroke)
	
	mtime_unixfile = -Inf
	# A variable controlling whether UNIX time information should be read from files in the event:
	readFiles = TRUE
	if(length(unixfile)>0){
		mtime_unixfile = file.info(unixfile)$mtime
		# Update if old function was used:
		if(mtime_unixfile < "2016-07-23 20:00:00 CEST"){
			warnings("Replacing old UNIX_time.TSD file with new and compressed file (all files prior to 2016-07-23 20:00:00 CEST are replaced)")
			fresh = TRUE
			}
		}
		
		
	### Get the UNIX_time data of the event in memory if present, saving time for large events: ###
	# Check for the UNIX_time information in memory, and use this if the time of the file in that information equals the time of the current UNIX_time file:
	# Accept if the difference in time is less than one second:
	getUNIX_timeInMemory = !isTRUE(fresh) && length(unixfile)>0 && saveData && exists(UNIX_time_mtime) && exists(UNIX_time_data) && abs(get(UNIX_time_mtime) - mtime_unixfile) < 1
	if(getUNIX_timeInMemory){
		out = get(UNIX_time_data)
		readFiles = FALSE
		file = FALSE
		if(identical(t,"all")){
			t = seq_along(out$f000)
			}
		}	
	else{
		# Using recursive=TRUE discards empty directories:
		event = list.files_caseInsensitive(event, full.names=TRUE, recursive=recursive)
		# Remove the files "UNIX_time.tsd" and "sgPM.tsd" from the list if present, as these only contain information about the unix time points of the other files and the segmentation parameter values, respectively:
		event = rmUnix_Sgpm(event)
		nfiles = length(event)
	
		if(nfiles==0){
			warning(paste("Event \"",event,"\" is empty of valid files (not UNIX_time.tsd and sgPM.tsd)",sep=""))
			return()
			}
	
		# Identify the tsd directory, if present:
		if(length(grep("tsd", tolower(tail(pathparts(dir), 1))))==0){
			tsddir = gregexpr("/tsd",dir)[[1]]
			if(tsddir[1]>0){
				attsd = max(tsddir)
				dir = substr(dir,1,attsd+4)
				}
			}
	
		# Treat all the files if t=="all":
		if(identical(t,"all")){
			t = seq_len(nfiles)
			}
	
		
		##### Execution and output #####
		out = list()
		# Get modification times:
		mtime_event = file.info(event)$mtime
		
		# Try reading the "UNIX_time.tsd" file:
		if(!isTRUE(fresh) && length(unixfile)>0){
			# Read the "UNIX_time.tsd" file:
			out <- read.TSD(unixfile, t="all", dimension=FALSE, header=FALSE, info=FALSE, drop.out=FALSE, use.raw=FALSE)
			
			# If information is given for only one file, re-list the data in 'out' (due to read.TSD, which drops dimensions if only one time step is requested):
			if(!is.list(out$l000)){
				out = lapply(out, list)
				}
			# Unlist f000 and t000 for convenience:
			out$f000 = unlist(out$f000, use.names=FALSE)
			out$t000 = unlist(out$t000, use.names=FALSE)
			
			out$I000 = out$I000[[1]]
			out$U000 = out$U000[[1]]
			out$t000 = sapply(out$t000, as.POSIXct, origin="1970-01-01 00:00.00 UTC")
			
			# Expand compactly stored i000:
			areseq = sapply(out$i000, function(xx) length(xx)==3 && head(xx, 1) == 0)
			out$i000[areseq] = lapply(out$i000[areseq], function(xx) seq.int(xx[2], xx[3]))
			if(length(out$I000)==3 && head(out$I000, 1)==0){
				out$I000 = seq.int(out$I000[2], out$I000[3])
				}
			# Add u000:
			out$u000 = lapply(out$i000, function(xx) out$U000[xx])
			
			if(saveData){
				assign(UNIX_time_data, out, envir=.GlobalEnv)
				assign(UNIX_time_mtime, mtime_unixfile, envir=.GlobalEnv)
				}
			
			# Check whether the filenames given in the unix time file are the same as the actual filenames, and whether the "UNIX_time.tsd" file is the latest:
			### names_unixfile = sort(tolower(c(out$f000)))[t]
			### names_filelist = sort(tolower(event))[t]
			names_unixfile = tolower(c(out$f000))[t]
			names_filelist = tolower(event)[t]
			### if(length(names_unixfile)==length(names_filelist)){
			### 	parts_names_unixfile1 = strsplit(names_unixfile[1],"/|\\\\")[[1]]
			### 	parts_names_filelist1 = strsplit(names_filelist[1],"/|\\\\")[[1]]
			### 	if(sum(parts_names_unixfile1!=parts_names_filelist1) < 2){
			### 		identicalnames = TRUE
			### 		}
			### 	else{
			### 		identicalnames = FALSE
			### 		}
			### 	}
			### else{
			### 	identicalnames = FALSE
			### 	}
			
			#identicalnames <- identical(names_unixfile==names_filelist)
			identicalnames <- identical(names_unixfile, names_filelist)
			latestfile = !any(difftime(mtime_unixfile, mtime_event[t])<0)
			
			# Determine whether a new UNIX time file should be written:
			if((allow.old && identicalnames) || (identicalnames && latestfile)){
				file = FALSE
				readFiles = FALSE
				}
			else if(length(file)>0 && !is.character(file[1])){
				file = NULL
				}
			}
		}
	
	# If fresh==NULL, return 'out' before trying to get the time information (used in noise.path()):
	if(length(fresh)==0){
		return(out)
		}
	
	###################################################################################################
	##### Get the new time information of file for which UNIX_time information has not been read: #####
	###################################################################################################
	if(!identical(readFiles,FALSE)){
		
		# Discard files that have not been changed:
		eventClean = setdiff(event, out$f000)
		nfilesClean = length(eventClean)
	
		# Read the UNIX_time information of the changed or new files:
		if(length(eventClean)>0){
			data = suppressWarnings(read.TSDs(eventClean, var=c("utim", "labl"), t="all", clean=FALSE, keep.all=TRUE, cores=cores))
			data = lapply(data, unlist, use.names=FALSE)
			utim_all = data[names(data)=="utim"]
			labl_all = data[names(data)=="labl"]
			rm(data)
			# Warning if label information is not present:
			if(identical(eventClean, event) && length(unlist(lapply(labl_all, length), use.names=FALSE))==0){
				warning("The event contains no TSD files")
				}
			
				# Discard time information for ctd-files:
				arectd = tolower(file_ext(eventClean))=="ctd"
				if(any(arectd)){
					utim_all[arectd] = vector("list", sum(arectd))
				}
				# Add the information that was read for the new or changed files:
				out$u000 = c(out$u000, utim_all)
				out$l000 = c(out$l000, labl_all)
				out$f000 = c(out$f000, eventClean)
			}
			
		
		
		
		# Remove information in the current UNIX_time file for files that are no longer present:
		present = out$f000 %in% event
		if(any(present)){
			out[validvarList] = lapply(out[validvarList], function(x) x[present])
		}
		
		
		
		
		# Merge with existing imfomation:
		
		### (1) File name f000:
		# Get the indexes 'i000' of the sorted unique unix time points:
		uniqueutim = sort(unique(unlist(out$u000, use.names=FALSE)))
		out$i000 = lapply(out$u000, function(x) as.double(findInterval(x, uniqueutim)))
		# Insert ones in the empty elements:
		NULLelments = which(unlist(lapply(out$i000,length), use.names=FALSE)==0)
		out$i000[NULLelments] = as.list(ones(length(NULLelments)))
		
		# Insert the unique 'utim' and 'indt' (U000, I000):
		out$U000 = list(sort(unique(unlist(out$u000, use.names=FALSE))))
		out$I000 = list(sort(unique(unlist(out$i000, use.names=FALSE))))
		
		# Then compress all sequences to c(head(x,1), NA, tail(x,1)) in i000 and I000:
		out$i000_full = out$i000
		areseq = sapply(out$i000, function(xx) diff(range(xx))==length(xx)-1 && length(xx) > 3)
		out$i000[areseq] = lapply(out$i000[areseq], function(xx) c(0, range(xx)))
		
		out$I000_full = out$I000
		if(diff(range(out$I000))==length(out$I000)-1 && length(out$I000) > 3){
			out$I000 = c(0, range(out$I000))
			}
		
		
		# Order by the order in 'event':
		orderEvent = match(event, out$f000)
		out[validvarList] = lapply(out[validvarList], function(x) x[orderEvent])
		
		# Insert modification time:
		out$t000 = as.character(mtime_event)
		}
	
	# Default is to write to a file named "UNIX_time.tsd" in the event:
	if(any(length(file)==0,nchar(file)==0,isTRUE(file))){
		file = file.path(dir,"UNIX_time.tsd")
		}	
	# Write the UNIX_time file:
	if(!identical(file,FALSE) && sum(unlist(lapply(out,length), use.names=FALSE))>0){
		listsTooDeep = lapply(out, function(x) if(length(x)>0) which(sapply(x,is.list)) else FALSE)
		listsTooDeep = unique(unlist(listsTooDeep, use.names=FALSE))
		if(length(listsTooDeep)>0){
			warning(paste("The following files are improperly organized (list at one or more of the time steps, which is not supported). Most likely the time information is not present for all time steps:\n",paste(event[listsTooDeep],collapse="\n")))
			}
		write.TSD(out[varToWrite], file, header=list(dtyp=processdtyp[varToWrite]), dimension=TRUE, numt=nfiles, keep.null=TRUE)
		}
	
	# Add the file names as names to the list elements for each variable:
	if("f000" %in% names(out)){
		addfilenames = c("i000", "l000", "t000")
		for(i in which(names(out) %in% addfilenames)){
			if(length(out[[i]])>0){
				names(out[[i]]) = out$f000
				}
			}
		}
	
	# Re-insert the full time step indices, and remove compressed:
	if(length(out)>7){
		out$i000 = out$i000_full
		out$I000 = out$I000_full
		out$i000_full = NULL
		out$I000_full = NULL
		}
	# Return the subset specified by 't':
	out[c("f000", "i000", "u000", "l000", "t000")] = lapply(out[c("f000", "i000", "u000", "l000", "t000")],function(x) x[t])
	# Unlist the I000 and U000 which were put in a list during writing to indicate one single time step for these variables:
	out[c("I000", "U000")] = lapply(out[c("I000", "U000")], unlist)
	
	
	invisible(out)
	##################################################
	##################################################
	}
