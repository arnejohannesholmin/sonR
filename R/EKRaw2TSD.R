#*********************************************
#*********************************************
#' Convert from Simrad or LSSS Profos files to TSD files.
#'
#' \code{EKRaw2TSD} Converts from raw data files in the EKRaw format to TSD files, readable by R using the functions read.TSD, resd.TSDs and read.event.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory. If given as a full path to the directory of the raw data, the value of 'esnm' is interpreted from the data. Can also be given as a vector of files directly.
#' @param filenr  is the indices of the files to convert.
#' @param t  is a vector of the time steps to process in each file.
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param eventname  is the name of the event, by default deduced from the file path.
#' @param CTD_station  is the path to the directory holding the CTD files. If not given (default) the directory is attemptedly set to file.path(Acoustics_datasets_directory(),"CTDs").
#' @param event_tsd  is the path to the TSD event, by default deduced from the file path.
#' @param endian  is the endian of the file, defaulted to .Platform$endian (changed from "big" by Arne Johannes Holmin 2012-07-31).
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param timeOffset  is the time offset of the datagram.
#' @param minTimeDiff  is the minimum difference in time betwee two pings.
#' @param drop  is TRUE to drop dimensions of the data.
#' @param msg  is TRUE to print a time bar during reading.
#' @param prenumt  is the number of time steps to reserve in the temporary file sto which data are saved during reading in order to restrict memeory usage. High values result in fewer temporary file but higher CPU time. The default value of 10 should be relatively optimal.
#' @param hpr0  is the pressure at sea level, which must be given externally. Defaulted to 101325 pascal.
#' @param Pain  is TRUE if pressure "ihpr" is given in Pascal and FALSE if given in decibar relative to surface pressure (10000 Pascal, giving values approximately equivalent to water depth).
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#' @param psze  is the vertical position of the acoustic instrument, given as a negative value in meters, defaulted to -7.
#' @param variableRange  should be set to TRUE for fishery sonars, 
#' @param gain  optinal gain values to apply to the data, typically from calibration.
#' @param sacr  optinal sacorrection values to apply to the data, typically from calibration.
#' @param correctTime  is TRUE to correct the data for offset in the time of the pings, determined by the difference in the time of the first ping and the first NMEA time information.
#' @param ow  is FALSE to not overwrite existing TSD files.
#' @param dira_offset  is an optional offset for the azimuth angle of the beams (useful to get the correct angle of fishery sonars).
#' @param write  is FALSE to only return the data and not write to TSD file.
#' @param cali  is FALSE to omit any calibration data located in the directory "Calibration" located at the top level of the cruise (alongside "Events"), or the path to a file if preferred.
#' @param toTS  is TRUE to apply the TS calibration instead for the Sv calibration.
#' @param na.rm  is TRUE to remove missing pings.
#' @param cleanNMEA  is FALSE to return the data as read, without shaving off incomplete time steps, 1 to remove incomplete and duplicate time steps, and 2 to additionally clean missing info at the end.
#' @param mergeFiles  is FALSE to not merge beams, vessel and rawvessel files at the end of the funciton, but rather keep the individual files.
#' @param cores  is an integer specifying the number of cores to run the compression over in parallel (should be lower than the number of cores in the computer).
#' @param tres  The time resolution of the compressed data in seconds.
#' @param xres  The sailed distance resolution of the compressed data in meters.
#' @param zres  The depth resolution of the compressed data in meters.
#' @param rres  The range resolution of the compressed data in meters.
#' @param bres  The beam resolution of the compressed data in integer number.
#' @param funvbsc  is the function to apply in the compression, either given as function or as a string, in which case the strings "mean" and "median" represents fast versions of the functions with the corresponding names (sum()/length() and fastMedian(), respectively).
#' @param funt  is the same as funvbsc, but used for averaging vessel data in the new time/distance bins.
#' @param adds  is a list of additional data overriding corresponding variables in 'data'
#' @param split used in psx.TSD().
#' @param skipAngles  is TRUE to discard electircal angles from the data (saves time).
#' @param origin  is either the time index of the origin, or the origin itself, given as c(longitude, latitude).
#' @param filesize  is the maximum size of the merged files.
#' @param chunksize  is the maximum size of the chunks of file read at the time.
#' @param apply.range.offset  Logical:: If TRUE strip the data of samples by the range offset Ro (for raw=1).
#' @param thr1m Logical: If TRUE apply a rule that TVG closer than 1 m should not cause increasedd level (as used in LSSS).
# 
# LSSS2TSD:
#
# @param event  is the identifier of the event,  given as the path to the event holding the folders "Export" and "Reports".
# @param filesize  is the maximum size of the merged files.
# @param chunksize  is the maximum size of the chunks of file read at the time.
#' @param keep.temp  is TRUE to delete the temporary files from individual school files.
# @param cores  is an integer specifying the number of cores to run the compression over in parallel (should be lower than the number of cores in the computer).
#' @param ...  For EKRaw2TSD(), arguments passed to combine.TSD() and compr.TSD(), and for LSSS2TSD(), inputs complementing the variables extracted from the files. Several of these are missing in the data as per 2016, such as (numt = number of time steps, and dimensions are given in square brackets):
#' \describe{
#'	\item{"numb"}{Number of beams [numt], interpreted from the file names}
#'	\item{"lenb"}{Lengths of the beams [numb, numt], set to the maximum range of the data}
#'	\item{"psze"}{Position of the transducer below sea surface, negative value [numt], defaulted to -7}
#'	\item{"diraSpan"}{The span of the sonar, given by the direction of the first and last beam in radians, where 0 is to the starboard side [2], defaulted to c(0, 2*pi)}
#'	\item{"asps"}{Average speed of sound in m/s [numt], iterated from the range values given the initial guess defaulted to 1480}
#'	\item{"bwtx"}{Beam width horizontally in radians [numb, numt], required}
#'	\item{"bwty"}{Beam width vertically in radians [numb, numt], required}
#'	\item{"bmmd"}{Beam modes, where 0: approximately horizontally oriented disc, like fishery sonar, 1: downward oriented echosounder, 2: downward oriented fan, 3: 3-D sonar, not vertical [numb, numt], required}
#'	\item{"freq"}{Frequency [numb, numt], sometimes required}
#'	\item{"eqba"}{Equivalent beam angle [numb, numt], not required}
#'	\item{"sacr"}{Sa-correction [numb, numt], not required}
#'	\item{"tpow"}{Transmit power [numt], not required}
#'	\item{"gai1", "gai1"}{Gain values [numt], not required}
#'	}
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom parallel makeCluster parLapply stopCluster
#' @importFrom SimradRaw readcalfile readEKRaw
#' @importFrom TSD combine.TSD read.TSD read.TSDs write.TSD
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @export
#' @rdname EKRaw2TSD
#'
EKRaw2TSD <- function(event, filenr="all", t="all", cruise=NULL, esnm="MS70", eventname=NULL, CTD_station=NULL, event_tsd=NULL, event_sadv=NULL, endian="little", dir.data=NULL, timeOffset=0, minTimeDiff=Inf, drop=TRUE, msg=TRUE, prenumt=10, hpr0=NULL, Pain=FALSE, TVG.exp=2, psze=-7, variableRange=TRUE, gain=NULL, sacr=NULL, correctTime=FALSE, ow=TRUE, dira_offset=0, write=TRUE, cali=TRUE, toTS=FALSE, na.rm=TRUE, cleanNMEA=1, mergeFiles=TRUE, cores=1, 
	# Inputs used in compr.TSD:
	tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("median","mean"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, 
	# Inputs used elsewhere:
	filesize=3e8, chunksize=1e8, apply.range.offset=FALSE, thr1m=FALSE, ...){
	
	############### LOG: ###############
	# Start: 2012-07-07 - Clean version.
	# Update: 2015-03-11 - Fixed bug with NMEA strings.

	
	##################################################
	##################################################
	########## Preparation ##########
	# Get file names for the last ping of the previous file group, the first ping of the current, and the rest of the pings:
	getFileNames123 <- function(x){
		file_path_sans_extx = file_path_sans_ext(x)
		file_extx <- file_ext(x)
		x1 <- file.path(paste0(file_path_sans_extx, "_1."), file_extx)
		x2 <- file.path(paste0(file_path_sans_extx, "_2."), file_extx)
		x3 <- file.path(paste0(file_path_sans_extx, "_3."), file_extx)
		list(x1, x2, x3, x)
	}
	
	compress <- any(sapply(c(tres, xres, zres, rres, bres), length)>0)
	if(isTRUE(write)){
		write <- c("p", "b", "v", "rv", "c")
	}
	
	# Check for the existence of the event of raw files:
	event_raw <- event.path(event=event, cruise=cruise, esnm=esnm, dir.data=dir.data, dir.type="raw")
	if(!file.exists(event_raw$event)){
		warning("Event not correctly specified")
		return()
	}
	if(length(eventname)==0){
		eventname <- event_raw$eventname
	}
	event_without_filetype <- dirname(event_raw$event)
	if(length(event_tsd)==0){
		event_tsd <- file.path(event_without_filetype, "tsd")
	}
	
	# List of raw files of the event. Accept a list of files, and not only a directory:
	if(identical(file.info(event[1])$isdir, FALSE)){
		filelist <- event
	}
	else{
		filelist <- list.files(event_raw$event,full.names=TRUE)
		filelist <- filelist[file_ext(filelist) == "raw"]
	}
	rawdir <- dirname(filelist[1])
	
	# If there are no raw files in the directory the function ends here:
	if(length(filelist) == 0){
		warning(paste("No raw files in the given directory ",event_raw$event,sep=""))
		return()
	}
	# File numbers:
	if(length(filenr)==0 || identical(filenr, "all")){
		filenr <- seq_along(filelist)
	}
	else if(!is.numeric(filenr)){
		warning("filenr must be numeric or \"all\"")
	}
	# Do not run more cores than pings, while keeping the possibly integer type of 'cores':
	if(length(filenr)<cores){
		cores <- filenr
	}
		
	# 'event_tsd' is name of the directory for the TSD files:
	if(!file.exists(event_tsd)){
		suppressWarnings(dir.create(event_tsd, recursive=TRUE))
	}
	
	
	
	########### Processing: ###########
	##### Moving through the files of the event: #####
	# Check if all the files are present in the tsd-directory, in which case the function is terminated:
	thisd <- readEKRaw(filelist[1], t=1:2, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, drop.out=drop, msg=FALSE, prenumt=prenumt)
	esnm <- EKRaw2TSD_getesnm(thisd, 1)
	
	# Create filenames:
	pingsfiles <- file.path(event_tsd, gsub(".raw",".pings",basename(filelist)))
	pingsfiles <- getFileNames123(pingsfiles)
	pingsfiles1 <- pingsfiles[[1]]
	pingsfiles2 <- pingsfiles[[2]]
	pingsfiles3 <- pingsfiles[[3]]
	pingsfiles <- pingsfiles[[4]]
	# Create beams-filenames:
	beamsfile <- file.path(event_tsd, paste(eventname,"_",esnm,".beams", sep=""))
	beamsfiles <- file.path(event_tsd, gsub(".raw",".beams",basename(filelist)))
	# Create vessel-filenames:
	vesselfile <- file.path(event_tsd, paste(eventname,"_",esnm,".vessel", sep=""))
	vesselfiles <- file.path(event_tsd, gsub(".raw",".vessel",basename(filelist)))
	# Create rawvessel-filenames:
	rawvesselfile <- file.path(event_tsd, paste(eventname,"_",esnm,"_rawVessel.tsd", sep=""))
	rawvesselfiles <- file.path(event_tsd, gsub(".raw","_rawVessel.tsd",basename(filelist)))
	# Create ctd-filenames:
	ctdfile <- file.path(event_tsd, paste(eventname,"_",esnm,".ctd", sep=""))
	ctdfiles <- file.path(event_tsd, gsub(".raw",".ctd",basename(filelist)))
	
	if(!ow){
		existingfiles <- list.files(event_tsd, recursive=TRUE, full.names=TRUE)
		if(length(existingfiles)){
			cat("Existing files: \n",paste0(existingfiles, collapse="\n"), "\n")
			cat("New files: \n",paste0(c(pingsfiles, beamsfile, vesselfile, rawvesselfile, ctdfile), collapse="\n"), "\n")
			ans <- readline("Overwrite existing files? (y/n)")
			if(!tolower(ans)=="y"){
				cat("Files already exist. Not overwriting\n")
				return()
			}
		}
	}
		
	# Look for a calibration file:
	if(is.character(cali) && file.exists(cali)){
		calfiles <- cali
		cali <- TRUE
	}
	else if(!identical(cali, FALSE)){
		# Look for calibration files or a folder called "calibration" holding calibration files:
		# Calibration files are xml files and must contain the string "cali" in the file name
		
		# Split the path to the directory holding the raw files into individual folder names:
		s <- strsplit(rawdir, .Platform$file.sep)[[1]]
		# Expand the paths to all the partent folders:
		ss <- sapply(seq_along(s), function(x) paste(s[seq_len(x)], collapse=.Platform$file.sep))
		# Get all directories, and search for those named "calibration" (case insinsitive). Order so that the directories closest to the rawdir comes first:
		dirs <- rev(unlist(lapply(ss, list.dirs, recursive=FALSE)))
		dirsbase <- basename(dirs)
		validDirs <- c(rawdir, dirs[grep("calibration", dirsbase, ignore.case=TRUE)])
		calfiles <- list.files(validDirs, full.names=TRUE)
		
		### calfiles <- list.files(file.path(dirname(dirname(dirname(dirname(dirname(filelist[1]))))), "Calibration"), full.names=TRUE, recursive=TRUE)
		### # Only read xml files:
		### calfiles <- calfiles[tolower(file_ext(calfiles))=="xml"]
		### calfiles <- calfiles[grep(esnm[1], calfiles)]
	}
	else{
		calfiles <- NULL
	}
	calfiles <- calfiles[tolower(file_ext(calfiles))=="xml"]
	# Old requirement of sonar name in the calibration file name:
	#calfiles <- calfiles[grep(esnm[1], basename(calfiles), ignore.case=TRUE)]
	
	# Add calibration data given in the input:
	if(!is.list(cali)){
		cali <- list()
	}
	# Read the calibration file:
	if(length(calfiles)>0){
		if(length(calfiles)>1){
			### if(is.character(cali)){
			### 	calfiles <- calfiles[grep(tolower(cali), tolower(calfiles))]
			### }
			### else{
			### 	warning(paste0("Several calibration files found for the system ", esnm[1], ". First chosen:", calfiles[1]))
			### 	calfiles <- calfiles[1]
			### }
			calfiles <- calfiles[1]
			#warning(cat("Multiple calibration files detected. The closest chosen", sep=""))
		}
		cat("Data calibrated from the following file (disregarding frequency):", calfiles, "\n")
		cali <- c(cali, readcalfile(calfiles))
	}
	else{
		if(is.sonar(esnm=esnm, fishery=TRUE)){
			cali <- c(cali, list(rofs=3))
			warning(paste0("Data of event ", rawdir, " not calibrated, maybe calibration file is missing for the system? The range offset 'Ro' due to the delay in signal processing set to the default:", cali$rofs))
		}
	}
			
	# Move through the list of raw files and read and possibly compress:
	if(cores>1 || identical(cores, 1L)){
		# Generate the clusters of time steps:
		cl<-makeCluster(cores)
		# Read all the .pings files, and merge at each time step:
		cat("Parallel processing on",cores,"cores:\n")
		fileinds <- parLapply(cl, filenr, EKRaw2TSD_oneFile_write, filelist=filelist, pingsfiles=pingsfiles, vesselfiles=vesselfiles, rawvesselfiles=rawvesselfiles, beamsfiles=beamsfiles, pingsfiles1=pingsfiles1, pingsfiles2=pingsfiles2, pingsfiles3=pingsfiles3, prenumt=prenumt, t=t, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, drop=drop, msg=msg, na.rm=na.rm, correctTime=correctTime, TVG.exp=TVG.exp, dira_offset=dira_offset, compress=compress, write=write, cali=cali, toTS=toTS, psze=psze, cleanNMEA=cleanNMEA, 
			# Inputs used in compr.TSD:
			tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, apply.range.offset=apply.range.offset, thr1m=thr1m, ...)
		# End the parallel processing:
		stopCluster(cl)
	}
	else{
		for(i in filenr){
			data <- EKRaw2TSD_oneFile_write(i=i, filelist=filelist, pingsfiles=pingsfiles, vesselfiles=vesselfiles, rawvesselfiles=rawvesselfiles, beamsfiles=beamsfiles, pingsfiles1=pingsfiles1, pingsfiles2=pingsfiles2, pingsfiles3=pingsfiles3, prenumt=prenumt, t=t, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, drop=drop, msg=msg, na.rm=na.rm, correctTime=correctTime, TVG.exp=TVG.exp, dira_offset=dira_offset, compress=compress, write=write, cali=cali, toTS=toTS, psze=psze, cleanNMEA=cleanNMEA, 
			# Inputs used in compr.TSD:
			tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, apply.range.offset=apply.range.offset, thr1m=thr1m, ...)
		}
	}
	
	
	if(compress || mergeFiles){
		# Read the vessel files, check for duplicated data, and merge:
		vessel <- read.TSDs(vesselfiles, t="all", clean=FALSE, merge=TRUE)
		vesselDuplicated <- duplicated(vessel$mtim)
	}
	else{
		vessel <- list()
	}
	# If there are any duplicated time steps, these are due to compression. Identify the duplicated time steps, and merge the compressed data by a weighted mean (not completely accurate in case median is used during compression):
	if(compress){
		if(any(vesselDuplicated)){
			# Merge the last ping of one file and the first ping of the next for all consecutive pairs of files:
			for(i in seq_len(length(pingsfiles1)-1)){
				pings3 <- read.TSD(pingsfiles3[i])
				pings1 <- read.TSD(pingsfiles1[i+1])
				if(pings3$mtim==pings1$mtim){
					for(j in seq_along(pings3)){
						pings3[[i]] <- (pings3[[i]] + pings1[[i]]) / 2
					}
					write.TSD(pingsfiles3[i], pings3, numt=1)
					unlink(pingsfiles1[i+1])
					pingsfiles1 <- pingsfiles1[-(i+1)]
				}
			}
		}
		
		# Simply merge the pings files
		allpingsfiles <- sort(c(pingsfiles1, pingsfiles2, pingsfiles3))
		### mergedfiles <- combine.TSD(allpingsfiles, reserve=TRUE, recursive=recursive, filesize=filesize, chunksize=chunksize, ...)$x_merged
		mergedfiles <- combine.TSD(allpingsfiles, reserve=TRUE, recursive=TRUE, filesize=filesize, chunksize=chunksize, ...)
		
		# The move the merged files to the original directory:
		newmergedfiles <- file.path(dirname(pingsfiles1[1]), basename(mergedfiles))
		
		file_path_sans_extx <- file_path_sans_ext(newmergedfiles)
		file_extx <- file_ext(newmergedfiles)
		newmergedfiles <- file.path(paste0(substr(file_path_sans_extx, nchar(newmergedfiles[,2])-2), "."), file_extx)
	
		#newmergedfiles <- get.ext(newmergedfiles, parts.out=TRUE)
		#newmergedfiles <- file.path(newmergedfiles[,1], paste0(substr(newmergedfiles[,2], 1, nchar(newmergedfiles[,2])-2), ".", newmergedfiles[,3]))
		if(any(duplicated(newmergedfiles))){
			warning("Duplicated new pings filenames, maybe due to too low value of 'filesize'")
		}
		file.rename(mergedfiles, newmergedfiles)
		# Remove the old files:
		unlink(allpingsfiles)
	}
	else{
		vessel$indt <- seq_along(vessel$mtim)
	}
	
	
	ctd <- NULL
	rawvessel <- NULL
	if(mergeFiles){
		# Read the beams files and merge:
		beams <- read.TSDs(beamsfiles, t="all", clean=FALSE, merge=TRUE)
		beams$indt <- seq_along(beams$mtim)
		# Read the rawvessel files and merge:
		rawvessel <- read.TSDs(rawvesselfiles, t="all", clean=FALSE, merge=TRUE)
		# Read the vessel files and merge:
		vessel <- read.TSDs(vesselfiles, t="all", clean=FALSE, merge=TRUE)
		vessel$indt <- seq_along(vessel$mtim)
		
		# Try adding missing log ('sadv') from the corresponding echosounder event if present:
		if((length(vessel$sadv) == 0 || all(is.na(vessel$sadv))) && length(event_sadv) && is.character(event_sadv) && file.exists(event_sadv)){
			# Look for the .vessel file:
			l <- list.files(event_sadv, "*.vessel")
			if(length(l)){
				vessel_sadv <- read.event(event=event_sadv, var="vessel", t="all")
				# Add the sadv from the echosounder event:
				vessel <- addVesselFromEchoosunder(vessel, echosounder=vessel_sadv)
			}
		}
	
	
		# Write the remaining data regardless of compression of the acoustic data:
		if("b" %in% write){ 
			write.TSD(beams, beamsfile)
			unlink(beamsfiles)
		}
		if("v" %in% write){ 
			write.TSD(vessel, vesselfile, header=list(dtyp=list(mtim="doub", lonv="doub", latv="doub", sadv="doub")))
			unlink(vesselfiles)
		}
		if("rv" %in% write){ 
			write.TSD(rawvessel, rawvesselfile, numt=1, header=list(dtyp=list(imtm="doub", ilnv="doub", iltv="doub", isdv="doub")))
			unlink(rawvesselfiles)
		}
	}
	if("c" %in% write){
		ctd <- CTD2TSD(CTD_station, outfile=ctdfile, hpr0=hpr0, Pain=Pain, vessel=vessel)
	}
	
	# Return individual lists or merge to one list:
	#invisible(c(data, vessel, rawvessel, ctd))
	
	return(event_tsd)
}
#'
#' @export
#' @rdname EKRaw2TSD
#'
EKRaw2TSDs <- function(event, esnm=c("SU90", "EK60"), dir.type=c("tsd", "Work"), cores=1){
	
	# Identify echosounder and run this first:
	isSBE <- sonR_implemented(esnm, "SBE")
	esnm <- c(esnm[isSBE], esnm[!isSBE])
	isSBE <- sonR_implemented(esnm, "SBE")
	isOFS <- sonR_implemented(esnm, "OFS")
	
	# Add the raw directory to the dir.type:
	dir.type <- c("raw", dir.type)
	
	# Get the full paths to the events:
	esnm_dir.type <- expand.grid(esnm=esnm, dir.type=dir.type)
	suppressWarnings(events <- apply(esnm_dir.type, 1, function(x) event.path(event, esnm=x["esnm"], dir.type=x["dir.type"])$event))
	events <- matrix(events, ncol=length(dir.type))
	dimnames(events) <- list(esnm, dir.type)
	events <- as.data.frame(events, stringsAsFactors=FALSE)
	
	# Run first the conversion from EKRaw to TSD, linking to the echosounder TSD directory, which should be the first directory as per the ordering of 'esnm' above:
	out <- lapply(events$raw, EKRaw2TSD, cores=cores, event_sadv=events$tsd[1])

	# Return a vector of paths to the TSD events:
	return(out)
	# Then convert the Work files:
	#lapply(events$Work, work2TSD, cores=cores)
}
#'
#' @export
#' @rdname EKRaw2TSD
#'
getSchoolsFromWork <- function(event, esnm="SU90", esnmLog=NULL, cores=1){
	
	# Use the sonar data to get the log per default:
	if(length(esnmLog) == 0){
		esnmLog <- esnm
	}
	
	# Get the full paths to the events:
	event_raw <- event.path(event, esnm=esnm, dir.type="raw")$event
	eventLog_tsd <- event.path(event, esnm=esnmLog, dir.type="tsd")$event
	
	# Get the school info and save a data frame with one row per school and vessel data appended to the columns:
	temp <- readLSSSWorkOFS(event=event_raw, cores=cores)

	# Get only the per school data:
	schools <- lapply(temp, "[[", "school")
	schools <- data.table::rbindlist(schools)
	schools <- as.data.frame(schools)

	# Add unix start and stop time:
	schools <- addDateTime(schools, prefix="Start")
	schools <- addDateTime(schools, prefix="Stop")

	# Add also mid Unix time:
	schools$utim <- (schools$Startutim + schools$Stoputim) / 2
	
	# Order by time:
	schools <- schools[order(schools$utim), ]

	##### Link the school data to the vessel data of the sonar (which has been linked with the echosounder data to get the log, which is named 'sadv' (sailed distance of the vessel) in TSD: #####
	# Read the vessel info from the echosounder:
	vesselSBE <- read.event(event=eventLog_tsd, var="vessel", t="all")

	#schools <- addVesselFromEchoosunder(schools, echosounder=vesselSBE, var=c("sadv", "lonv", "latv"))
	schools <- addVesselFromEchoosunder(schools, echosounder=vesselSBE)
	
	return(schools)
}
#'
#' @import data.table
#' @export
#' @rdname EKRaw2TSD
#'
aggregateSchoolsFromWork <- function(x, by="sadv", delta=1){
	x <- data.table::as.data.table(x)
	
	x$key <- floor(x[,by, with=FALSE] / delta)
	
	x$HeadingRad <- x$Heading * pi/180
	
	out <- x[,  .(
		"StartDateTime" = StartDateTime[1], 
		"StopDateTime" = StopDateTime[1], 
		"AverageSpeed" = mean(Speed, na.rm=TRUE), 
		"AverageHeading" = atan2(sum(Speed * sin(HeadingRad), na.rm=TRUE), sum(Speed * cos(HeadingRad), na.rm=TRUE)), 
		"SumArea" = sum(Area, na.rm=TRUE), 
		"AverageArea" = mean(Area, na.rm=TRUE), 
		"AverageAlongBeamSize" = mean(AlongBeamSize, na.rm=TRUE), 
		"AverageAlongRingSize" = mean(AlongRingSize, na.rm=TRUE), 
		"NumSchools" = length(unique(Id)), 
		"AverageDepth" = mean(Depth, na.rm=TRUE), 
		"MedianSv" = median(Mean, na.rm=TRUE)
		), by="key"]
		
	data.table::setnames(out, "key", by)
		
	out <- as.data.frame(out)	
	
	return(out)
}

#'
#' @importFrom TSD labl.TSD write.TSD
#' 
EKRaw2TSD_oneFile_write <- function(i, filelist, pingsfiles, vesselfiles, rawvesselfiles, beamsfiles, pingsfiles1=NULL, pingsfiles2=NULL, pingsfiles3=NULL, prenumt=10, t="all", endian="little", timeOffset=0, minTimeDiff=Inf, msg=TRUE, na.rm=TRUE, correctTime=FALSE, TVG.exp=2, dira_offset=0, compress=FALSE, write=TRUE, cali=TRUE, toTS=FALSE, psze=-7, cleanNMEA=1, 
	# Inputs used in compr.TSD:
	tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("median","mean"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, apply.range.offset=FALSE, thr1m=FALSE, ...){

	# Declare the variable names:
	beamsnames <- TSD::labl.TSD("EKRaw2TSD_b")
	vesselnames <- TSD::labl.TSD("EKRaw2TSD_v")
	rawvesselnames <- TSD::labl.TSD("EKRaw2TSD_r")
	pingsnames <- TSD::labl.TSD("EKRaw2TSD_p")

	# Read the data from the raw file:
	data <- EKRaw2TSD_oneFile(i=i, filelist=filelist,  prenumt=prenumt, t=t, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, msg=msg, na.rm=na.rm, correctTime=correctTime, TVG.exp=TVG.exp, dira_offset=dira_offset, cali=cali, toTS=toTS, psze=psze, skipAngles=skipAngles, cleanNMEA=cleanNMEA, apply.range.offset=apply.range.offset, thr1m=thr1m)
	
	if(length(data)==0){
		rm(data)
		gc()
		return(i)
	}
	numt <- length(data$mtim)
		
	# Compress the data:
	if(compress){
		data <- compr.TSD(data, tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, ...)
	}
		
	# Write the first time step, the last time step and the center time steps, so that the last and first time steps can be merged between data coming from consecutive raw files:
	if("p" %in% write){ 
		if(compress){
			# Write the first ping:
			TSD::write.TSD(data[pingsnames], pingsfiles1[i], t=1, numt=numt, header=list(dtyp=list(vbsc="floa")))
			# Write the middle pings:
			TSD::write.TSD(data[pingsnames], pingsfiles2[i], t=seq(2, numt-1), numt=numt, header=list(dtyp=list(vbsc="floa")))
			# Write the last ping:
			TSD::write.TSD(data[pingsnames], pingsfiles3[i], t=numt, numt=numt, header=list(dtyp=list(vbsc="floa")))
		}
		else{
			TSD::write.TSD(data[pingsnames], pingsfiles[i], numt=numt, header=list(dtyp=list(vbsc="floa")))
		}
	}
	
	# Write the remaining data regardless of compression of the acoustic data:
	if("b" %in% write){ 
		TSD::write.TSD(data[beamsnames], beamsfiles[i], numt=numt)
	}
	if("v" %in% write){ 
		TSD::write.TSD(data[vesselnames], vesselfiles[i], numt=numt, header=list(dtyp=list(mtim="doub", lonv="doub", latv="doub", sadv="doub")))
	}
	if("rv" %in% write){ 
		TSD::write.TSD(data[rawvesselnames], rawvesselfiles[i], numt=1, header=list(dtyp=list(imtm="doub", ilnv="doub", iltv="doub", isdv="doub")))
	}
	
	rm(data)
	gc()
	#invisible(data)
	i
}
#'
#' @importFrom SimradRaw NMEA2vessel readEKRaw readEKRaw_power2sv.TSD readEKRaw_stripNA getRangeOffsetInUnitsOfSamples
#' @importFrom TSD listOfEqual2array mergeListKeepDimensions NAs
#' @importFrom utils tail head
#' @importFrom SimradRaw offsetSamples
#'
EKRaw2TSD_oneFile <- function(i, filelist,  prenumt=10, t="all", endian="little", timeOffset=0, minTimeDiff=Inf, msg=TRUE, na.rm=TRUE, correctTime=FALSE, TVG.exp=2, dira_offset=0, cali=TRUE, toTS=FALSE, psze=-7, skipAngles=TRUE, cleanNMEA=1, apply.range.offset=FALSE, thr1m=FALSE){

	# Function for treating the beam directions and beam widths:
	EKRwa2TSD_treat_beams <- function(beams, dira_offset){
		# Convert beam directions and widths to the standard definitions used in the TSD format:
		beams <- Simrad_dir(beams, dira_offset=dira_offset)
		beams <- Simrad_bwt(beams)
		beams
	}
	# Function for treating the vessel dynamics:
	EKRwa2TSD_treat_vessel <- function(vessel, rawvessel){
		# Discard time steps which does not follow chronologically:
		imtmDiff <- diff(rawvessel$imtm)
		rawvessel <- rawvessel[unlist(lapply(rawvessel,length))>0]
		rawvessel$irzv <- (rawvessel$irzv * pi/180) %% (2*pi)
		# Get vessel dynamics from rawvessel
		vessel <- lapply(vessel, mergeListKeepDimensions)
		vessel <- vessel[unlist(lapply(vessel,length))>0]
		
		# Interpolate the raw vessel data to the ping times:
		#vessel[c("rtzv", "latv", "lonv", "ispv", "sadv")] <- lapply(rawvessel[c("irzv", "iltv", "ilnv", "iisv", "isdv")], function(xx) if(length(xx)) Hmisc::approxExtrap(rawvessel$imtm, xx, vessel$mtim)$y else NAs(length(vessel$mtim)))
		vessel[c("rtzv", "latv", "lonv", "ispv", "sadv")] <- lapply(rawvessel[c("irzv", "iltv", "ilnv", "iisv", "isdv")], function(xx) if(length(xx)==1) rep(xx, length.out=length(vessel$mtim)) else if(length(xx)>1) approx(rawvessel$imtm, xx, vessel$mtim, rule=2)$y else NAs(length(vessel$mtim)))
		
		## Get positions of the time points in the data:
		#imtmIntervals <- c(rawvessel$imtm[1]-imtmDiff[1]/2, rawvessel$imtm[-length(rawvessel$imtm)] + imtmDiff/2, rawvessel$imtm[length(rawvessel$imtm)]+imtmDiff[length(imtmDiff)]/2)
		#atpingstime <- findInterval(vessel$mtim, imtmIntervals)
		#atpingstime[atpingstime==0] <- 1
		#
		## Assign the vessel data:
		#vessel[c("rtzv", "latv", "lonv", "ispv", "sadv")] <- lapply(rawvessel[c("irzv", "iltv", "ilnv", "iisv", "isdv")], function(xx) if(length(xx)) xx[atpingstime] else zeros(length(vessel$mtim)))
		#vessel$rtzv <- rawvessel$irzv[atpingstime]
		#vessel$latv <- rawvessel$iltv[atpingstime]
		#vessel$lonv <- rawvessel$ilnv[atpingstime]
		#vessel$ispv <- rawvessel$iisv[atpingstime]
		#vessel$sadv <- rawvessel$isdv[atpingstime]
		vessel
	}
	
	# Declare lists (except rawvessel, which is directly deduced from the NMEA string):
	pings <- list()
	beams <- list()
	vessel <- list()
	rawvessel <- list()
	ctd <- list()

	########## (1) Reading raw data file number 'f': ##########
	cat("Reading raw file number ",i," of ",length(filelist), " (",filelist[i],") ...\n",sep="")
	thisd <- try(suppressWarnings(readEKRaw(filelist[i], t=t, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, drop.out=FALSE, msg=msg, prenumt=prenumt, na.rm=na.rm)), silent=TRUE)
	
	if(length(thisd)==1 && is(thisd)=="try-error"){
		warning(paste0("Error when reading the Simrad raw file using readEKRaw(). Empty data returned: \"", thisd, "\""))
		return(NULL)
	}
	# If nothing was read, return empty lists:
	if(length(thisd)==0 || sum(unlist(lapply(thisd$data$pings, length)))==0){
		return(NULL)
	}
	
	# Correct the ping time with the offset between the first NMEA time and the first ping time:
	timediff=NA
	if(correctTime){
		NMEAimtm <- NMEA2vessel(thisd$data$NMEA$string, cleanNMEA=cleanNMEA)
		if(length(NMEAimtm$imtm)){
			# This addition was found ad hoc for the MS70 data using the script "Time error in MS70 .R":
			if(strff("ms70", EKRaw2TSD_getesnm(thisd, numt=1))){
				TOV_offset <- -0.25/86400
		}
			else{
				TOV_offset <- 0
		}
			timediff <- NMEA2vessel(thisd$data$NMEA$string, cleanNMEA=cleanNMEA)$imtm[1] - thisd$data$pings$time[1] + TOV_offset
			thisd$data$pings$time <- thisd$data$pings$time + timediff
		}
	}
	
	## If specified, split the channels into separate time steps with beam modes. This is only applied if 'bmmd' is given ???????? When is this ever used ??????? (2016-12-12):
	#if(length(bmmd) == thisd$header$transceivercount){
	#	head2 <- function(x, n=1, N=2){
	#		if(NCOL(x)==N){
	#			head(t(x), n)
	#		}
	#		else{
	#			head(x, n)
	#		}
	#	}
	#	numb <- thisd$header$transceivercount
	#	thisd <- lapply(thisd$data, function(x1) lapply(x1, function(x2) if(is.list(x2)) lapply(x2, head2, 1, numb) else head2(x2, 1, numb)))
	#}
	
	########## (2) Store beams variables: ##########
	numt <- NCOL(thisd$data$pings$time)
	# (2.1) esnm - the name of the system:
	beams$esnm <- EKRaw2TSD_getesnm(thisd, numt)
	# (2.3) mtim - Matlab time:
	beams$mtim <- thisd$data$pings$time[1,]
	# (2.4) asps - average speed of sound:
	beams$asps <- thisd$data$pings$soundvelocity[1,]
	# (2.5) numb - number of beams:
	numb <- if(length(thisd$header$transceivercount)) thisd$header$transceivercount else length(thisd$data$config$frequency)
	beams$numb <- rep(numb, numt)
	# (2.6) indi - beam indices:
	beams$indi <- matrix(seq_len(beams$numb[1]), beams$numb[1], numt)
	# (2.7) freq - frequency:
	beams$freq <- thisd$data$pings$frequency
	# (2.8) absr - absorption coefficient:
	beams$absr <- thisd$data$pings$absorptioncoefficient
	# (2.9) sint - sampling interval duration
	beams$sint <- thisd$data$pings$sampleinterval[1,]
	# (2.10) rres - radial resolution
	beams$rres <- beams$asps * beams$sint / 2
	# (2.11) plsl - sound pulse duration:
	beams$plsl <- thisd$data$pings$pulselength[1,]
	# (2.12) plsl - sound pulse duration:
	if(length(thisd$data$pings$transducerdepth) && all(thisd$data$pings$transducerdepth != 0, na.rm=TRUE)){
		beams$psze <- -thisd$data$pings$transducerdepth[1,]
	}
	else{
		beams$psze <- rep(psze, numt)
	}
	# (2.12) lenb - lengths of the beams:
	beams$lenb <- thisd$data$pings$count

	# Raw1:
	if(length(thisd$data$pings$beamwidthhorizontaltx)){
		# (2.13) dirx - direction (x) of the beams:
		beams$dirx <- thisd$data$pings$dirx
		# (2.14) diry - direction (y) of the beams:
		beams$diry <- thisd$data$pings$diry
		# (2.15) dirz - direction (z) of the beams:
		beams$dirz <- thisd$data$pings$dirz
		# (2.18) bwtl - beam width along ship:
		beams$bwtl <- thisd$data$pings$beamwidthalongshiprx
		# (2.19) bwtt - beam width athwart ship:
		beams$bwtt <- thisd$data$pings$beamwidthathwartshiprx
		# (2.20) bwth - beam width horizontal (used in SX90):
		beams$bwth <- thisd$data$pings$beamwidthhorizontaltx
		# (2.21) bwtv - beam width vertical (used in SX90:
		beams$bwtv <- thisd$data$pings$beamwidthverticaltx
		# (2.24) eqba - equivalent beam angle:
		beams$eqba <- thisd$data$pings$equivalentbeamangle
		# (2.25) sacr - sa-correction:
		beams$sacr <- thisd$data$pings$sacorrection
		# (2.26) tpow - transmit power:
		beams$tpow <- thisd$data$pings$transmitpower
		# (2.27) gai1 - gain at emission (used in SX90):
		beams$gai1 <- thisd$data$pings$gaintx
		# (2.28) gai2 - gain at reception (used in SX90):
		beams$gai2 <- thisd$data$pings$gainrx
		# (2.29) gain - gain:
		beams$gain <- beams$gaintx + beams$gainrx
	}
	# Raw0:
	else{
		# Get the position in the pulselengthtable:
		atPulselength <- which(thisd$data$config$pulselengthtable == beams$plsl[1], arr.ind=TRUE)
		# (2.13) dirx - direction (x) of the beams:
		beams$dirx <- matrix(thisd$data$config$dirx, beams$numb[1], numt)
		# (2.14) diry - direction (y) of the beams:
		beams$diry <- matrix(thisd$data$config$diry, beams$numb[1], numt)
		# (2.15) dirz - direction (z) of the beams:
		beams$dirz <- matrix(thisd$data$config$dirz, beams$numb[1], numt)
		# (2.18) bwtl - beam width along ship:
		beams$bwtl <- matrix(thisd$data$config$beamwidthalongship, beams$numb[1], numt)
		# (2.19) bwtt - beam width athwart ship:
		beams$bwtt <- matrix(thisd$data$config$beamwidthathwartship, beams$numb[1], numt)
		# (2.24) eqba - equivalent beam angle:
		beams$eqba <- matrix(thisd$data$config$equivalentbeamangle, beams$numb[1], numt)
		# (2.25) sacr - sa-correction:
		#beams$sacr <- matrix(thisd$data$config$sacorrectiontable[,1], beams$numb[1], numt)
		beams$sacr <- matrix(thisd$data$config$sacorrectiontable[atPulselength], beams$numb[1], numt)
		# (2.26) tpow - transmit power:
		beams$tpow <- thisd$data$pings$transmitpower
		# (2.29) gain - gain:
		#beams$gain <- matrix(thisd$data$config$gain, beams$numb[1], numt)
		beams$gain <- matrix(thisd$data$config$gaintable[atPulselength], beams$numb[1], numt)
	}
	
	# (2.30) bmmd - beam mode:
	if(length(thisd$data$pings$beammode)>0){
		beams$bmmd <- thisd$data$pings$beammode[1,]
		bmmdFromDira <- apply(beams$dirx, 2, function(xx) any(abs(xx-xx[1])>1e-3)) * 2
		if(any(beams$bmmd != bmmdFromDira)){
			warning(paste("Some beam modes (bmmd) did not agree with the beam angles, and were changed:\n",paste(which(beams$bmmd != bmmdFromDira), collapse=", ")))
			beams$bmmd <- bmmdFromDira
		}
	}
	
						
	########## (2) Store raw vessel variables: ##########
	NMEAstring <- thisd$data$NMEA$string
	rawvessel <- NMEA2vessel(NMEAstring, cleanNMEA=cleanNMEA)
	
	
	########## (3) Store per ping vessel variables: ##########
	# (3.2) mtim - Matlab time:
	vessel$mtim <- thisd$data$pings$time[1,]
	# (3.3) pszv - Heave:
	vessel$pszv <- thisd$data$pings$heave[1,]
	# (3.4) rtxv - Pitch:
	vessel$rtxv <- (thisd$data$pings$pitch[1,] * pi/180) %% (2*pi)
	# (3.5) rtyv - Roll:
	vessel$rtyv <- (thisd$data$pings$roll[1,] * pi/180) %% (2*pi)
	# (3.10) terr - Ping time error:
	vessel$terr <- rep(timediff*86400, numt)
	
	
	########## (4) Store pings variables: ##########
	# (3.2) mtim - Matlab time:
	pings$mtim <- thisd$data$pings$time[1,]
	# (3.4) angl - Electrical angle along ship:
	if(!skipAngles && length(thisd$data$pings$alongship_e)){
		pings$angl <- thisd$data$pings$alongship_e
		# No longer used, since the new readEKRaw() outputs an array [#samples, #beams, #pings] padded with NAs.
		#pings$angl <- TSD::listOfEqual2array(pings$angl)
	}
	# (3.5) angt - Electrical angle athwart ship:
	if(!skipAngles && length(thisd$data$pings$alongship_e)){
		pings$angt <- thisd$data$pings$alongship_e
		# No longer used, since the new readEKRaw() outputs an array [#samples, #beams, #pings] padded with NAs.
		#pings$angt <- TSD::listOfEqual2array(pings$angt)
	}
	
	# (3.3) vbsc - Volume backscattering coefficient:
	isOmniSonar <- tolower(beams$esnm[[1]]) %in% c("sx80", "sh80", "su80", "sx90", "sh90", "su90")
	temp <- readEKRaw_power2sv.TSD(thisd$data$pings$power, beams, cali=cali, tiltcorr=isOmniSonar, toTS=toTS, list.out=TRUE)
	pings$vbsc <- temp$vbsc
	beams$Cgai <- temp$Cgai
	beams$Csac <- temp$Csac
	beams$Ctcr <- temp$Ctcr
	beams$Ccal <- temp$Ccal
	# Also add the range offset, either from the calibration (in the raw1 case) or by the raw0 value:
	beams$rofs <- rep(if(isOmniSonar) cali$rofs else beams$asps * beams$plsl / 4, length.out=numt)
	# Get the sample offset due to the range offset:
	beams$sofs <- rep(getRangeOffsetInUnitsOfSamples(beams), length.out=numt)
	
	if(apply.range.offset){
		# Strip the data of samples by the range offset Ro (for raw=1):
		pings$vbsc <- offsetSamples(pings$vbsc, beams=beams)
		# Add TVG:
		pings$vbsc <- apply.TVG(pings$vbsc, beams=beams, linear=TRUE, TVG.exp=TVG.exp, Ro=0, thr1m=thr1m)
	}
	else{
		# Add TVG:
		pings$vbsc <- apply.TVG(pings$vbsc, beams=beams, linear=TRUE, TVG.exp=TVG.exp, thr1m=thr1m)
	}
	
	
	rm(thisd)
	rm(temp)
	gc()
	# Add lengths of beams and number of beams:
	pings$lenb <- beams$lenb
	pings$numb <- beams$numb
	
	########## (5) Write pings file: ##########
	# Remove missing pings:
	if(na.rm && any(is.na(pings$mtim))){
		valid <- is.na(pings$mtim)
		beams <- lapply(beams, readEKRaw_stripNA, valid)
		vessel <- lapply(vessel, readEKRaw_stripNA, valid)
		pings <- lapply(pings, readEKRaw_stripNA, valid)
	}
	
	beams <- EKRwa2TSD_treat_beams(beams, dira_offset)
	
	if(length(rawvessel$imtm)){
		vessel <- EKRwa2TSD_treat_vessel(vessel, rawvessel)
	}
	
	# Output the data read from the Simrad raw file:
	c(pings, beams, vessel, rawvessel)
}

EKRaw2TSD_getesnm <- function(thisd, numt){
	
	# Due to a change in the header, 'soundername' can also be named 'ApplicationName':
	if(!length(thisd$header$soundername) && length(thisd$header$ApplicationName)){
		esnm <- thisd$header$ApplicationName
	}
	else if(length(thisd$header$soundername)){
		esnm <- thisd$header$soundername
	}
	else{
		stop("Neither 'soundername' nor 'ApplicationName' included in the header:\n\t", paste(names(thisd$header), thisd$header, collapse=", ", sep=" = "))
	}
	
	if( length(grep("MBS", esnm) > 0) || length(grep("MS70", esnm) > 0) || identical(thisd$header$transceivercount, 500) ){
		rep("MS70", numt)
	}
	else if( length(grep("MBES", esnm) > 0) || length(grep("ME70", esnm) > 0) ){
		rep("ME70", numt)
	}
	else if( length(grep("ER60", esnm) > 0) || length(grep("EK60", esnm) > 0) ){
		rep("EK60", numt)
	}
	else if( length(grep("SX80", esnm) > 0) ){
		rep("SX80", numt)
	}
	else if( length(grep("SH80", esnm) > 0) ){
		rep("SH80", numt)
	}
	else if( length(grep("SU80", esnm) > 0) ){
		rep("SU80", numt)
	}
	else if( length(grep("SX90", esnm) > 0) ){
		rep("SX90", numt)
	}
	else if( length(grep("SH90", esnm) > 0) ){
		rep("SH90", numt)
	}
	else if( length(grep("SU90", esnm) > 0) ){
		rep("SU90", numt)
	}
	else{
		rep(esnm, numt)
	}
}

addVesselFromEchoosunder <- function(x, echosounder, var="sadv"){
	
	# The UNIX time must be present in both lists:
	x$utim <- utim.TSD(x)
	echosounder$utim <- utim.TSD(echosounder)
	if(!length(x$utim) || !length(echosounder$utim)){
		stop("Both lists must contain time coercable to UNIX time using utim.TSD() utim")
	}
	
	indt <- findInterval(x$utim, echosounder$utim)
	# If there are times in x$utim before the first time in echosounder$utim, replace the zeros in 'indt' with 1 (assign these to the first echosounder interval) and add a warning:
	are0 <- which(indt == 0)
	if(length(are0)){
		indt[are0] <- 1
		warning("Some time steps of the sonar Work files are before the first time of the echosounder. These time steps are assigned to the first time of the echosounder when adding log form the echosounder to the data from the Work files. These are the following ping numbers in the sonar Work files: ", prettyIntegers(are0))
	}
	
	# Add log:
	for(this in var){
		x[[this]] <- echosounder[[this]][indt]
	}
	#x$sadv <- echosounder$sadv[indt]
	
	return(x)
}




addDateTime <- function(x, prefix="", tz="UTC"){
	DateVar <- paste0(prefix, "Date")
	TimeVar <- paste0(prefix, "Time")
	
	temp <- list()
	temp$DateTime <- as.POSIXct(paste(x[[DateVar]], x[[TimeVar]]), tz=tz)
	temp$utim <- unclass(temp$DateTime)
	
	names(temp) <- paste0(prefix, names(temp))
	
	x <- cbind(x, temp)
	
	return(x)
}


# Function for reading work files:
readLSSSWorkOFS_one <- function(x){
	
	# Function for a character variable of a list, splitting it by the 'split' argument, and returning the resulting numeric vector:
	getVarNumeric <- function(x, name, split=NULL){
		out <- x[[name]]
		if(length(split)){
			out <- strsplit(out, split)[[1]]
		}
		out <- as.numeric(out)
		out
	}
	
	# The segmentation mask is stored by a start range followed the number of voxels included in the mask, then optionally th enumber of voxels excluded and then the number incldued, and so on:
	expandBeams <- function(ping){
		
		# Function for expanding a segmentation mask og one vector from the compact form (a starting index, then the number of incldued indices, followed by pairs of number of excluded and included indices) to a vector of indices of the mask:
		expandSamples <- function(x){
			# The data are saved by a starting index, then the number of incldued indices, followed by pairs of number of excluded and included indices:
			starts <- cumsum(x)
			ends <- starts[-1] - 1
			starts <- starts[-length(starts)]
			# Pick out only the included:
			included <- TSD::odd(seq_along(starts))
			starts <- starts[included]
			ends <- ends[included]
			
			mask <- unlist(lapply(seq_along(starts), function(i) seq.int(starts[i], ends[i])))
			mask
		}
		
		# Function for reading and expanding one beam of the segmentation mask:
		expandOneBeam <- function(oneBeam){
			beamInd <- getVarNumeric(oneBeam, "channel")
			sampleInd <- getVarNumeric(oneBeam, "range", split=" ")
			sampleInd <- expandSamples(sampleInd)
			out <- cbind(sampleInd, beamInd)
			out
		}
		
		# Get the number of samples along the beams, the number of beams, and the actual segmentation mask:
		lenb <- getVarNumeric(ping$.attrs, "samples")
		numb <- getVarNumeric(ping$.attrs, "beams")
		beams <- ping[names(ping) == "beam"]
		
		# Expand all the beams and rbind and convert to flat indices:
		arr.ind <- lapply(beams, expandOneBeam)
		arr.ind <- do.call(rbind, arr.ind)
		ind <- TSD::arr.ind2ind(arr.ind, c(lenb, numb))
		ind
	}
	
	# Read the xml file (possibly from a local file):
	out <- Rstox::downloadXML(x)
	
	
	# Get info per ping:
	getInfo <- function(x){
		as.list(unlist(unname(x)))
	}
	
	# Get the info per ping:
	info <- t(sapply(out$aspects$perPing, getInfo))
	dimNamesInfo <- dimnames(info)
	dimInfo <- dim(info)
	# Unlist and convert to a matrix:
	info <- unlist(info)
	info <- array(info, dim=dimInfo)
	dimnames(info) <- dimNamesInfo
	
	# Convert to data frame:
	numericCols <- !is.na(suppressWarnings(sapply(head(info, 1), as.numeric)))
	info <- as.data.frame(info, stringsAsFactors=FALSE)
	info[,numericCols] <- apply(info[,numericCols], 2, as.numeric)
	
	
	
	
	#info$DateTime <- as.POSIXct(paste(info$Date, info$Time), tz="UTC")
	#info$utim <- unclass(info$DateTime)
	
	info <- addDateTime(info)
	
	info$indtInFile <- as.numeric(basename(info$file))
	info$fileID <- dirname(info$file)
	
	
	school <- getInfo(out$aspects$AllPings)
	#names(school) <- paste("School", names(school), sep="_")
	
	# The Work files have the string "N/A" for NAs, so we need to convert to NA before trying to assign to numeric:
	convertToNA <- function(x, NAstring="N/A"){
		replace(x, grep(NAstring, x), NA)
	}
	school <- lapply(school, convertToNA)
	
	areNumeric <- function(x){
		!any(is.na(as.numeric(x)) && !is.na(x))
	}
	
	numeric <- suppressWarnings(sapply(school, areNumeric))
	school[numeric] <- lapply(school[numeric], as.numeric)
	
	# Get the school IDs:
	#schoolID <- as.numeric(out$aspects$AllPings$Id$Id)
	
	# Get the time and convert to UNIX time:
	#Time <- lapply(out$aspects$perPing, "[[", "Time")
	#Time <- sapply(Time, paste, collapse=" ")
	#utim <- TSD::ftim2utim(Time)
	
	## Get file ID and ping index in the file:
	#fileInd <- unlist(lapply(out$aspects$perPing, "[[", "File"))
	#indtInFile <- as.numeric(basename(fileInd))
	#fileID <- dirname(fileInd)
	
	# Combine all time and file info in a data frame:
	#time <- data.frame(
	#	Time = Time, 
	#	utim = utim, 
	#	fileID = fileID, 
	#	indtInFile = indtInFile, 
	#	stringsAsFactors = FALSE
	#)
	
	# Read the segmentation mask:
	mask <- lapply(out$mask, expandBeams)
	
	# Return the data:
	out <- list(
		Id = info$Id[1], 
		WorkFile = basename(x), 
		info = info, 
		school = school, 
		mask = mask
	)
	
	return(out)
}

readLSSSWorkOFS <- function(event, filenr="all", esnm="SU90", cores=1){
	# Locate work files:
	event_work <- sonR::event.path(event, dir.type="Work", esnm=esnm)$event
	l <- list.files(event_work, full.names=TRUE)
	
	# Read only the file with file extension "work":
	l <- l[tolower(tools::file_ext(l)) %in% "work"]
	
	if(is.numeric(filenr)){
		l <- l[filenr]
	}
	
	system.time(sv <- TSD::papply(l, readLSSSWorkOFS_one, cores=cores))
	invisible(sv)
}




# Function for converting work files to TSD files of segmentation masks:
work2TSD <- function(event, cores=1){
	
	# FINISH THIS LATER
	
	event_tsd <- sonR::event.path(event, dir.type="tsd")$event
	
	# Read all work files:
	sv <- readLSSSWorkOFS(event=event, cores=cores)
	
	# Merge segmentation masks for each ping, adding school ID as sgID:
	sgsc <- lapply(sv, "[[", "mask")
	numt <- sapply(sgsc, length)
	sgsc <- unlist(sgsc, recursive=FALSE)
	numx <- lengths(sgsc)
	sgID <- sapply(sv, "[[", "schoolID")
	sgID <- rep(sgID, numt)
	sgID <- lapply(seq_along(sgID), function(i) rep(sgID[i], numx[i]))
	
	
	utim <- unlist(lapply(sv, "[[", c("time", "utim")), recursive=FALSE)
	file <- unlist(lapply(sv, "[[", c("time", "fileID")), recursive=FALSE)
	pind <- unlist(lapply(sv, "[[", c("time", "indtInFile")), recursive=FALSE)
	
	# Split by time
	sgsc <- split(sgsc, utim)
	sgID <- split(sgID, utim)
	file <- split(file, utim)
	pind <- split(pind, utim)
	utim <- sort(unique(utim))
	
	# Collapse for each ping:
	sgscFlat <- lapply(sgsc, unlist, use.names=FALSE, recursive=FALSE)
	sgIDFlat <- lapply(sgID, unlist, use.names=FALSE, recursive=FALSE)
	
	
	# Read the unix time of the TSD files:
	u <- UNIX_time(event_tsd)
	
	# Select only the horizontal beams:
	bmmd <- read.event(event_tsd, var="bmmd", t="all")$bmmd
	atHorizontal <- which(bmmd == 0)
	
	
	indt <- atHorizontal[findInterval(utim, u$U000[atHorizontal])]
	
	# Do not accept duplicated indices:
	validPings <- !duplicated(indt)
	
	# Write the valid pings:
	seglist <- list(
		sgsc = sgscFlat[validPings], 
		sgID = sgIDFlat[validPings], 
		utim = u$U000[indt[validPings]]
	)
	
	segfile <- file.path(event_tsd, "profosSegmentationVertical.seg")
	write.TSD(seglist, segfile)
	
	
}







