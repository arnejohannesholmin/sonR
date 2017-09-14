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
#' @param correctTime  is FALSE to avoid correcting the data for offset in the time of the pings, determined by the difference in the time of the first ping and the first NMEA time information.
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
# 
# LSSS2TSD:
#
# @param event  is the identifier of the event,  given as the path to the event holding the folders "Export" and "Reports".
# @param filesize  is the maximum size of the merged files.
# @param chunksize  is the maximum size of the chunks of file read at the time.
#' @param keep.temp  is TRUE to delete the temporary files from individual school files.
# @param cores  is an integer specifying the number of cores to run the compression over in parallel (should be lower than the number of cores in the computer).
#' @param ...  For EKRaw2TSD(), arguments passed to merge_TSD() and compr.TSD(), and for LSSS2TSD(), inputs complementing the variables extracted from the files. Several of these are missing in the data as per 2016, such as (numt = number of time steps, and dimensions are given in square brackets):
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
#' @importFrom TSD merge_TSD read.TSD read.TSDs write.TSD
#' @importFrom tools file_ext file_path_sans_ext
#'
#' @export
#' @rdname EKRaw2TSD
#'
EKRaw2TSD<-function(event, filenr="all", t="all", cruise=NULL, esnm="MS70", eventname=NULL, CTD_station=NULL, event_tsd=NULL, endian="little", dir.data=NULL, timeOffset=0, minTimeDiff=Inf, drop=TRUE, msg=TRUE, prenumt=10, hpr0=NULL, Pain=FALSE, TVG.exp=2, psze=-7, variableRange=TRUE, gain=NULL, sacr=NULL, correctTime=TRUE, ow=TRUE, dira_offset=0, write=TRUE, bmmd=NULL, cali=TRUE, toTS=FALSE, na.rm=TRUE, cleanNMEA=1, mergeFiles=TRUE, cores=1, 
	# Inputs used in compr.TSD:
	tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("median","mean"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, 
	# Inputs used elsewhere:
	filesize=3e8, chunksize=1e8, apply.range.offset=FALSE, thr1m=FALSE, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-07 - Clean version.
	# Update: 2015-03-11 - Fixed bug with NMEA strings.

	
	##################################################
	##################################################
	########## Preparation ##########
	# Get file names for the last ping of the previous file group, the first ping of the current, and the rest of the pings:
	getFileNames123 <- function(x){
		file_path_sans_extx = file_path_sans_ext(x)
		file_extx = file_ext(x)
		x1 = file.path(paste0(file_path_sans_extx, "_1."), file_extx)
		x2 = file.path(paste0(file_path_sans_extx, "_2."), file_extx)
		x3 = file.path(paste0(file_path_sans_extx, "_3."), file_extx)
		list(x1, x2, x3, x)
		}
	compress = any(sapply(c(tres, xres, zres, rres, bres), length)>0)
	if(isTRUE(write)){
		write = c("p", "b", "v", "rv", "c")
		}
	
	# Check for the existence of the event of raw files:
	event_raw = event.path(event=event, cruise=cruise, esnm=esnm, dir.data=dir.data, dir.type="raw")
	if(!file.exists(event_raw$event)){
		warning("Event not correctly specified")
		return()
		}
	if(length(eventname)==0){
		eventname = event_raw$eventname
		}
	event_without_filetype = dirname(event_raw$event)
	if(length(event_tsd)==0){
		event_tsd = file.path(event_without_filetype, "tsd")
		}
	
	# List of raw files of the event. Accept a list of files, and not only a directory:
	if(identical(file.info(event[1])$isdir, FALSE)){
		filelist = event
		}
	else{
		filelist = list.files(event_raw$event,full.names=TRUE)
		filelist = filelist[file_ext(filelist) == "raw"]
		}
	rawdir <- dirname(filelist[1])
	
	# If there are no raw files in the directory the function ends here:
	if(length(filelist) == 0){
		warning(paste("No raw files in the given directory ",event_raw$event,sep=""))
		return()
		}
	# File numbers:
	if(length(filenr)==0 || identical(filenr, "all")){
		filenr = seq_along(filelist)
		}
	else if(!is.numeric(filenr)){
		warning("filenr must be numeric or \"all\"")
		}
	# Do not run more cores than pings, while keeping the possibly integer type of 'cores':
	if(length(filenr)<cores){
		cores = filenr
		}
		
	# 'event_tsd' is name of the directory for the TSD files:
	if(!file.exists(event_tsd)){
		suppressWarnings(dir.create(event_tsd, recursive=TRUE))
		}
	
	
	
	########### Processing: ###########
	##### Moving through the files of the event: #####
	# Check if all the files are present in the tsd-directory, in which case the function is terminated:
	thisd <- readEKRaw(filelist[1], t=1:2, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, drop.out=drop, msg=FALSE, prenumt=prenumt)
	esnm = EKRaw2TSD_getesnm(thisd, 1)
	
	# Create filenames:
	pingsfiles = file.path(event_tsd, gsub(".raw",".pings",basename(filelist)))
	pingsfiles = getFileNames123(pingsfiles)
	pingsfiles1 = pingsfiles[[1]]
	pingsfiles2 = pingsfiles[[2]]
	pingsfiles3 = pingsfiles[[3]]
	pingsfiles = pingsfiles[[4]]
	# Create beams-filenames:
	beamsfile = file.path(event_tsd, paste(eventname,"_",esnm,".beams", sep=""))
	beamsfiles = file.path(event_tsd, gsub(".raw",".beams",basename(filelist)))
	# Create vessel-filenames:
	vesselfile = file.path(event_tsd, paste(eventname,"_",esnm,".vessel", sep=""))
	vesselfiles = file.path(event_tsd, gsub(".raw",".vessel",basename(filelist)))
	# Create rawvessel-filenames:
	rawvesselfile = file.path(event_tsd, paste(eventname,"_",esnm,"_rawVessel.tsd", sep=""))
	rawvesselfiles = file.path(event_tsd, gsub(".raw","_rawVessel.tsd",basename(filelist)))
	# Create ctd-filenames:
	ctdfile = file.path(event_tsd, paste(eventname,"_",esnm,".ctd", sep=""))
	ctdfiles = file.path(event_tsd, gsub(".raw",".ctd",basename(filelist)))
	
	if(!ow){
		existingfiles = list.files(event_tsd, recursive=TRUE, full.names=TRUE)
		if(length(existingfiles)){
			cat("Existing files: \n",paste0(existingfiles, collapse="\n"), "\n")
			cat("New files: \n",paste0(c(pingsfiles, beamsfile, vesselfile, rawvesselfile, ctdfile), collapse="\n"), "\n")
			ans = readline("Overwrite existing files? (y/n)")
			if(!tolower(ans)=="y"){
				cat("Files already exist. Not overwriting\n")
				return()
				}
			}
		}
		
	# Look for a calibration file:
	if(is.character(cali) && file.exists(cali)){
		calfiles = cali
		cali = TRUE
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
		
		### calfiles = list.files(file.path(dirname(dirname(dirname(dirname(dirname(filelist[1]))))), "Calibration"), full.names=TRUE, recursive=TRUE)
		### # Only read xml files:
		### calfiles = calfiles[tolower(file_ext(calfiles))=="xml"]
		### calfiles = calfiles[grep(esnm[1], calfiles)]
		}
	else{
		calfiles = NULL
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
			### 	calfiles = calfiles[grep(tolower(cali), tolower(calfiles))]
			### }
			### else{
			### 	warning(paste0("Several calibration files found for the system ", esnm[1], ". First chosen:", calfiles[1]))
			### 	calfiles = calfiles[1]
			### }
			calfiles <- calfiles[1]
			#warning(cat("Multiple calibration files detected. The closest chosen", sep=""))
		}
		cat("Data calibrated from the following file (disregarding frequency):", calfiles, "\n")
		cali = c(cali, readcalfile(calfiles))
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
		fileinds = parLapply(cl, filenr, EKRaw2TSD_oneFile_write, filelist=filelist, pingsfiles=pingsfiles, vesselfiles=vesselfiles, rawvesselfiles=rawvesselfiles, beamsfiles=beamsfiles, pingsfiles1=pingsfiles1, pingsfiles2=pingsfiles2, pingsfiles3=pingsfiles3, prenumt=prenumt, t=t, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, drop=drop, msg=msg, na.rm=na.rm, correctTime=correctTime, bmmd=bmmd, TVG.exp=TVG.exp, dira_offset=dira_offset, compress=compress, write=write, cali=cali, toTS=toTS, psze=psze, cleanNMEA=cleanNMEA, 
			# Inputs used in compr.TSD:
			tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, apply.range.offset=apply.range.offset, thr1m=thr1m, ...)
		# End the parallel processing:
		stopCluster(cl)
		}
	else{
		for(i in filenr){
			data = EKRaw2TSD_oneFile_write(i=i, filelist=filelist, pingsfiles=pingsfiles, vesselfiles=vesselfiles, rawvesselfiles=rawvesselfiles, beamsfiles=beamsfiles, pingsfiles1=pingsfiles1, pingsfiles2=pingsfiles2, pingsfiles3=pingsfiles3, prenumt=prenumt, t=t, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, drop=drop, msg=msg, na.rm=na.rm, correctTime=correctTime, bmmd=bmmd, TVG.exp=TVG.exp, dira_offset=dira_offset, compress=compress, write=write, cali=cali, toTS=toTS, psze=psze, cleanNMEA=cleanNMEA, 
			# Inputs used in compr.TSD:
			tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, apply.range.offset=apply.range.offset, thr1m=thr1m, ...)
			}
		}
	
	
	if(compress || mergeFiles){
		# Read the vessel files, check for duplicated data, and merge:
		vessel = read.TSDs(vesselfiles, t="all", clean=FALSE, merge=TRUE)
		vesselDuplicated = duplicated(vessel$mtim)
		}
	else{
		vessel = list()
		}
	# If there are any duplicated time steps, these are due to compression. Identify the duplicated time steps, and merge the compressed data by a weighted mean (not completely accurate in case median is used during compression):
	if(compress){
		if(any(vesselDuplicated)){
			# Merge the last ping of one file and the first ping of the next for all consecutive pairs of files:
			for(i in seq_len(length(pingsfiles1)-1)){
				pings3 = read.TSD(pingsfiles3[i])
				pings1 = read.TSD(pingsfiles1[i+1])
				if(pings3$mtim==pings1$mtim){
					for(j in seq_along(pings3)){
						pings3[[i]] = (pings3[[i]] + pings1[[i]]) / 2
						}
						write.TSD(pingsfiles3[i], pings3, numt=1)
						unlink(pingsfiles1[i+1])
						pingsfiles1 = pingsfiles1[-(i+1)]
					}
				}
			}
		
		# Simply merge the pings files
		allpingsfiles = sort(c(pingsfiles1, pingsfiles2, pingsfiles3))
		### mergedfiles = merge_TSD(allpingsfiles, reserve=TRUE, recursive=recursive, filesize=filesize, chunksize=chunksize, ...)$x_merged
		mergedfiles = merge_TSD(allpingsfiles, reserve=TRUE, recursive=TRUE, filesize=filesize, chunksize=chunksize, ...)
		
		# The move the merged files to the original directory:
		newmergedfiles = file.path(dirname(pingsfiles1[1]), basename(mergedfiles))
		
		file_path_sans_extx = file_path_sans_ext(newmergedfiles)
		file_extx = file_ext(newmergedfiles)
		newmergedfiles = file.path(paste0(substr(file_path_sans_extx, nchar(newmergedfiles[,2])-2), "."), file_extx)
	
		#newmergedfiles = get.ext(newmergedfiles, parts.out=TRUE)
		#newmergedfiles = file.path(newmergedfiles[,1], paste0(substr(newmergedfiles[,2], 1, nchar(newmergedfiles[,2])-2), ".", newmergedfiles[,3]))
		if(any(duplicated(newmergedfiles))){
			warning("Duplicated new pings filenames, maybe due to too low value of 'filesize'")
			}
		file.rename(mergedfiles, newmergedfiles)
		# Remove the old files:
		unlink(allpingsfiles)
		}
	else{
		vessel$indt = seq_along(vessel$mtim)
		}
	
	
	ctd = NULL
	rawvessel = NULL
	if(mergeFiles){
		# Read the beams files and merge:
		beams = read.TSDs(beamsfiles, t="all", clean=FALSE, merge=TRUE)
		beams$indt = seq_along(beams$mtim)
		# Read the rawvessel files and merge:
		rawvessel = read.TSDs(rawvesselfiles, t="all", clean=FALSE, merge=TRUE)
		# Read the vessel files and merge:
		vessel = read.TSDs(vesselfiles, t="all", clean=FALSE, merge=TRUE)
		vessel$indt = seq_along(vessel$mtim)
	
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
		ctd = CTD2TSD(CTD_station, outfile=ctdfile, hpr0=hpr0, Pain=Pain, vessel=vessel)
		}
	
	# Return individual lists or merge to one list:
	invisible(c(data, vessel, rawvessel, ctd))
	##################################################
	##################################################
	}
