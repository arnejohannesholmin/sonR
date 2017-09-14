#*********************************************
#*********************************************
#' This funciton compresses the data stored in 'x' by averaging (either mean, trimmed mean or median) in depth and time bins. This is designed primarily for echosounder beams.
#'
#' @param data  is the list of inputs variables as returned from EKRaw2TSD().
#' @param tres  The time resolution of the compressed data in seconds.
#' @param xres  The sailed distance resolution of the compressed data in meters.
#' @param zres  The depth resolution of the compressed data in meters.
#' @param rres  The range resolution of the compressed data in meters.
#' @param bres  The beam resolution of the compressed data in integer number.
#' @param cores  is an integer specifying the number of cores to run the compression over in parallel (should be lower than the number of cores in the computer).
#' @param funvbsc  is the function to apply in the compression, either given as function or as a string, in which case the strings "mean" and "median" represents fast versions of the functions with the corresponding names (sum()/length() and fastMedian(), respectively).
#' @param funt  is the same as funvbsc, but used for averaging vessel data in the new time/distance bins.
#' @param adds  is a list of additional data overriding corresponding variables in 'data'
#' @param split used in psx.TSD().
#' @param skipAngles  is TRUE to discard electircal angles from the data (saves time).
#' @param origin  is either the time index of the origin, or the origin itself, given as c(longitude, latitude).
#' @param write  is FALSE to only return the data and not write to TSD file.
#' @param ...  further arguments passed to psx.TSD().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom parallel makeCluster parLapply stopCluster clusterSplit
#' @importFrom TSD labl.TSD merge_TSD read.TSD read.TSDs write.TSD
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils tail head
#' @importFrom stats weighted.mean
#'
#' @export
#' @rdname compr.TSD
#'
compr.event <- function(event, filenr="all", tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, cores=1, funvbsc=c("median","mean"), funt=c("median","mean"), adds=NULL, split=TRUE, skipAngles=TRUE, origin=1, z0=0, cruise=NULL, esnm="EK60", event_compr=NULL, ow=TRUE, msg=TRUE, write=TRUE, filesize=3e8, chunksize=3e8, keepEmpty=TRUE, clear_individual=TRUE, clear_along=FALSE, maxlenb=NULL, ...){
		
	########## Preparation ##########
	# Used in merge_TSD_by_Time():
	head1 = function(x, ...){
		head(x, 1, ...)
		}
	# Get file names for the last ping of the previous file group, the first ping of the current, and the rest of the pings:
	getComprFileNames123 <- function(fileBaseNamesSansExt, event_compr, ext="pings"){
		comprFileNames = file.path(event_compr, fileBaseNamesSansExt)
		comprFileNames1 = paste0(comprFileNames, "_1.", ext)
		comprFileNames2 = paste0(comprFileNames, "_2.", ext)
		comprFileNames3 = paste0(comprFileNames, "_3.", ext)
		list(
			comprFileNames1=comprFileNames1, 
			comprFileNames2=comprFileNames2, 
			comprFileNames3=comprFileNames3, 
			comprFileNames=comprFileNames)
	}
	
	compress = any(sapply(c(tres, xres, zres, rres, bres), length)>0)
	
	event = event.path(event=event, cruise=cruise, esnm=esnm)
	eventname = event$eventname
	
	if(length(event_compr)==0){
		# Get existing compression directories:
		comprDir = list.files(dirname(event$event), full.names=TRUE)
		comprDir = comprDir[grep("compr", basename(comprDir))]
		comprNr = suppressWarnings(as.numeric(gsub("(^.+compr_+)(\\d+)(_.+$)", "\\2", comprDir)))
		#comprNr = comprNr[!is.na(comprNr)]
		resAdd = list(tres, xres, zres, rres, bres)
		names(resAdd) = c("tres", "xres", "zres", "rres", "bres")
		resAdd = resAdd[sapply(resAdd, length)>0]
		resAdd = paste(names(resAdd), resAdd, sep="_", collapse="_")
		event_compr_final = file.path(dirname(event$event), paste0("tsd_compr_", max(comprNr, 0, na.rm=TRUE)+1, "_", resAdd))
		event_compr = paste0(event_compr_final, "_", "individual")
		}
	else{
		event_compr_final = NULL
	}
	
	# List of raw files of the event. Accept a list of files, and not only a directory:
	TIME = UNIX_time(event$event)
	filelist = unlist(TIME$f000, use.names=FALSE)
	# Get file base names of the pings files:
	arePingsFiles = which(file_ext(filelist) == "pings")
	indt = TIME$i000[arePingsFiles]
	#TIME[c("f000", "u000", "i000")] = lapply(TIME[c("f000", "u000", "i000")], "[", arePingsFiles)
	# Extract file names to be used for the compressed files:
	fileBaseNamesSansExt = file_path_sans_ext(basename(unlist(TIME$f000, use.names=FALSE)[arePingsFiles]))
	
	# If there are no raw files in the directory the function ends here:
	if(length(filelist) == 0){
		warning(paste("No pings files in the given directory ", event$event, sep=""))
		return()
		}
	# File numbers:
	if(length(filenr)==0 || identical(filenr, "all")){
		filenr = seq_along(fileBaseNamesSansExt)
		}
	else if(!is.numeric(filenr)){
		warning("filenr must be numeric or \"all\"")
		}
	# Do not run more cores than pings, while keeping the possibly integer type of 'cores':
	if(length(filenr)<cores){
		warning("The number of cores reduced to the number of pings-files")
		cores = length(filenr)
		}

	# Order 'filenr' chronologically:
	indtStart = lapply(indt[filenr], head, 1)
	lindtStart1 = sapply(indtStart, length)==1
	if(!all(lindtStart1)){
		warning("Some files did not have any time steps, and were excluded from the compression")
		filenr = filenr[lindtStart1]
		indtStart = lapply(indt[filenr], head, 1)
		}
	filenr = filenr[order(unlist(indtStart, use.names=FALSE))]

	# 'event_compr' is name of the directory for the TSD files (create it here and not before in case there are no pings files):
	if(!file.exists(event_compr)){
		suppressWarnings(dir.create(event_compr, recursive=TRUE))
		}
	
	# Create pings-filenames:
	pingsfiles = getComprFileNames123(fileBaseNamesSansExt, event_compr, ext="pings")
	# Create beams-filenames:
	beamsfiles = getComprFileNames123(fileBaseNamesSansExt, event_compr, ext="beams")
	# Create vessel-filenames:
	vesselfiles = getComprFileNames123(fileBaseNamesSansExt, event_compr, ext="vessel")
	
	
	########### Processing: ###########
	##### Moving through the files of the event: #####
	# Declare the variable names:
	beamsnames = labl.TSD("EKRaw2TSD_b")
	vesselnames = labl.TSD("EKRaw2TSD_v")
	pingsnames = labl.TSD("EKRaw2TSD_p")
	if(skipAngles){
		pingsnames = setdiff(pingsnames, c("angl", "angt"))
		}
	
	# Check if all the files are present in the tsd-directory, in which case the function is terminated:
	if(!ow){
		existingfiles = list.files(event, recursive=TRUE, full.names=TRUE)
		if(length(existingfiles)){
			#cat("Existing files: \n",paste0(existingfiles, collapse="\n"), "\n")
			#cat("New files: \n",paste0(c(pingsfiles[[4]], beamsfile, vesselfile, ctdfile), collapse="\n"), "\n")
			ans = readline("Overwrite existing files? (y/n)")
			if(!tolower(ans)=="y"){
				cat("Files already exist. Not overwriting\n")
				return()
				}
			}
		}
	
	# Move through the list of pings files and read and possibly compress:
	if(cores>1 || identical(cores, 1L)){
		# Generate the clusters of time steps:
		cores = min(cores, length(filenr))
		cl<-makeCluster(cores, outfile="")
		# Reorder the files so that they are read as closely in position on the drive as possible:
		cc = clusterSplit(cl, filenr)
		corenr = rep(seq_len(cores), sapply(cc, length))
		l = unlist(lapply(seq_len(cores), function(xx) seq(xx, by=cores, l=length(cc[[xx]]))), use.names=FALSE)
		corenr[rank(l)] = corenr
													################# Here we are creating a vector of files which is only as long as the length of filenr, whereas the indt is for all files, which may be unwise. The function will only work for filenr starting from 1 !!!!!!!!!!!!!!!!!!!!!!!!
		dumpfiles = NAs(length(indt))
		dumpfiles[filenr] = file.path(dirname(event_compr), paste0(basename(event_compr), "_core_", corenr, ".txt"))
		l = filenr[rank(l)]
		# Read all the .pings files, and merge at each time step:
		cat("Parallel processing on", cores, "cores:\n")
		# The output 'indices' is not used:
		indices = parLapply(cl, l, compr.event_oneFile_write, indt=indt, filelist=filelist, pingsfiles=pingsfiles, vesselfiles=vesselfiles, beamsfiles=beamsfiles, t="all", compress=compress, TIME=TIME, write=write, 
			# Inputs used in compr.TSD:
			tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, z0=z0, pingsnames=pingsnames, vesselnames=vesselnames, beamsnames=beamsnames, dumpfiles=dumpfiles, keepEmpty=keepEmpty, maxlenb=maxlenb, ...)
		# End the parallel processing:
		stopCluster(cl)
		}
	else{
		# The output 'indices' is not used:
		indices = NAs(length(filenr))
		for(i in filenr){
			indices[i] = compr.event_oneFile_write(i=i, indt=indt, filelist=filelist, pingsfiles=pingsfiles, vesselfiles=vesselfiles, beamsfiles=beamsfiles, t="all", compress=compress, TIME=TIME, write=write, 
			# Inputs used in compr.TSD:
			tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, z0=z0, pingsnames=pingsnames, vesselnames=vesselnames, beamsnames=beamsnames, keepEmpty=keepEmpty, maxlenb=maxlenb, ...)
			}
		}
	
	# If there are files in the final merged directory, move the last file to the temporarily-merged directory (but only if it has 1 time step), and apply merging by time steps:
	x_final = list.files(event_compr_final, full.names=TRUE)
	if(length(x_final)){
		# Split by file extension:
		ext_final = file_ext(x_final)
		# Get the final pings, beams, and vessel files:
		x_final = split(x_final, ext_final)
		names(x_final) = sort(unique(ext_final))
		x_final = x_final[names(x_final) %in% c("pings", "beams", "vessel")]
		# Check if the last files all have 1 time step, and move the files if so:
		lastfiles = sapply(x_final, tail, 1)
		lastNumt = read.TSDs(lastfiles, var="numt", clean=FALSE)
		lastNumt = lastNumt[names(lastNumt)=="numt"]
		if(all(sapply(lastNumt, length)==1)){
			file.rename(lastfiles, file.path(dirname(x_final[1]), basename(lastfiles)))
			}
		}
	
	# Merge the last file of the current file number and the first file of the next:
	pingsfiles = merge_TSD_by_Time(unlist(lapply(pingsfiles[1:3], "[", filenr), use.names=FALSE), fun=weighted.mean)
	beamsfiles = merge_TSD_by_Time(unlist(lapply(beamsfiles[1:3], "[", filenr), use.names=FALSE), fun=head1)
	vesselfiles = merge_TSD_by_Time(unlist(lapply(vesselfiles[1:3], "[", filenr), use.names=FALSE), fun=head1)
	#indt = seq_len(sum(read.TSDs(pingsfiles, var="numt", merge=TRUE, cores=cores)$numt))
	
	# Then merge the files using the linked option:
	#merge_TSD(pingsfiles, dir=event_compr_final, linked=list(beamsfiles, vesselfiles), adds=list(indt=indt), clear_along=TRUE)
	merge_TSD(pingsfiles, dir=event_compr_final, linked=list(beamsfiles, vesselfiles), clear_along=clear_along, skipLast=TRUE)
	
	# Rename all files by removing "_1" just before the file extension:
	l = list.files(event_compr_final, full.names=TRUE)
	file_path_sans_ext_l = file_path_sans_ext(l)
	newl = paste(substr(file_path_sans_ext_l, 1, nchar(file_path_sans_ext_l)-2), file_ext(l), sep=".")
	file.rename(l, newl)
	
	# Delete the individual files and rename the merged directory to the compressed directory:
	if(clear_individual){
		cat("Deleting individual files\n")
		unlink(event_compr, force=TRUE, recursive=TRUE)
		}
	
	# Return:
	invisible(list(event=event_compr_final, pingsfiles=pingsfiles, beamsfiles=beamsfiles, vesselfiles=vesselfiles))
	##################################################
	##################################################
	}
	
