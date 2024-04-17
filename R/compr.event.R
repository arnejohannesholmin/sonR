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
#' @param cores  is an integer specifying the number of cores to run the compression over in parallel (should be ler than the number of cores in the computer).
#' @param funvbsc  is the function to apply in the compression, either given as function or as a string, in which case the strings "mean" and "median" represents fast versions of the functions with the corresponding names (sum()/length() and fastMedian(), respectively). Default is mean, which is recommended for volume backscattering coefficient (sv) data, which are by nature exponentially distributed (backscattering from multiple targets), and using median will underestmate the true backscatter.
#' @param funt  is the same as funvbsc, but used for averaging vessel data in the new time/distance bins.
#' @param adds  is a list of additional data overriding corresponding variables in 'data'
#' @param split used in psx.TSD().
#' @param skipAngles  is TRUE to discard electircal angles from the data (saves time).
#' @param origin  is either the time index of the origin, or the origin itself, given as c(longitude, latitude).
#' @param write  is FALSE to only return the data and not write to TSD file.
#' @param ow  Logical: If TRUE, overwrite the latest compression (directory).
#' @param ...  further arguments passed to psx.TSD().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom parallel makeCluster parLapply stopCluster clusterSplit
#' @importFrom TSD labl.TSD combine.TSD read.TSD read.TSDs write.TSD
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils tail head
#' @importFrom stats weighted.mean
#'
#' @export
#' @rdname compr.event
#'
compr.event <- function(event, filenr="all", tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, cores=1, funvbsc=c("mean","median"), funt=c("median","mean"), adds=NULL, split=TRUE, skipAngles=TRUE, origin=1, z0=0, cruise=NULL, esnm="EK60", msg=TRUE, write=TRUE, filesize=3e8, chunksize=3e8, clear_individual=TRUE, clear_along=FALSE, maxlenb=NULL, ow=FALSE, ...){
	#compr.event <- function(event, filenr="all", tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, cores=1, funvbsc=c("mean","median"), funt=c("median","mean"), adds=NULL, split=TRUE, skipAngles=TRUE, origin=1, z0=0, cruise=NULL, esnm="EK60", event_compr=NULL, msg=TRUE, write=TRUE, filesize=3e8, chunksize=3e8, clear_individual=TRUE, clear_along=FALSE, maxlenb=NULL, ow=FALSE, ...){
		
	########## Preparation ##########
	# Used in merge_TSD_by_Time():
	head1 <- function(x, ...){
		head(x, 1, ...)
	}
	# Get file names for the last ping of the previous file group, the first ping of the current, and the rest of the pings:
	getComprFileNames123 <- function(fileBaseNamesSansExt, event_compr, ext="pings"){
		comprFileNames <- file.path(event_compr, fileBaseNamesSansExt)
		comprFileNames1 <- paste0(comprFileNames, "_1.", ext)
		comprFileNames2 <- paste0(comprFileNames, "_2.", ext)
		comprFileNames3 <- paste0(comprFileNames, "_3.", ext)
		list(
			comprFileNames1=comprFileNames1, 
			comprFileNames2=comprFileNames2, 
			comprFileNames3=comprFileNames3, 
			comprFileNames=comprFileNames)
	}
	
	compress <- any(sapply(c(tres, xres, zres, rres, bres), length)>0)
	
	event <- event.path(event=event, cruise=cruise, esnm=esnm)
	eventname <- event$eventname
	
	#if(length(event_compr)==0){
		# Get existing compression directories:
		comprDir <- list.files(dirname(event$event), full.names=TRUE)
		comprDir <- comprDir[grep("compr", basename(comprDir))]
		comprNr <- suppressWarnings(as.numeric(gsub("(^.+compr_+)(\\d+)(_.+$)", "\\2", comprDir)))
		
		if(length(comprNr) == 0){
			comprNr <- 1
		}
		else{
			comprNr <- max(comprNr, na.rm=TRUE) + if(!ow) 1 else 0
		}
		
		#comprNr <- comprNr[!is.na(comprNr)]
		resAdd <- list(tres, xres, zres, rres, bres)
		names(resAdd) <- c("tres", "xres", "zres", "rres", "bres")
		resAdd <- resAdd[sapply(resAdd, length)>0]
		resAdd <- paste(names(resAdd), resAdd, sep="_", collapse="_")
		event_compr_final <- file.path(dirname(event$event), paste0("tsd_compr_", comprNr, "_", resAdd))
		event_compr <- paste0(event_compr_final, "_", "individual")
		if(ow){
			if(file.exists(event_compr_final)){
				unlink(event_compr_final, force=TRUE, recursive=TRUE)
			}
			if(file.exists(event_compr)){
				unlink(event_compr, force=TRUE, recursive=TRUE)
			}
			
		}
	#}
	#else{
	#	event_compr_final <- NULL
	#}
	
	# List of raw files of the event. Accept a list of files, and not only a directory:
	TIME <- UNIX_time(event$event)
	filelist <- unlist(TIME$f000, use.names=FALSE)
	# Get file base names of the pings files:
	arePingsFiles <- which(file_ext(filelist) == "pings")
	indt <- TIME$i000[arePingsFiles]
	#TIME[c("f000", "u000", "i000")] = lapply(TIME[c("f000", "u000", "i000")], "[", arePingsFiles)
	# Extract file names to be used for the compressed files:
	fileBaseNamesSansExt <- file_path_sans_ext(basename(unlist(TIME$f000, use.names=FALSE)[arePingsFiles]))
	
	# If there are no raw files in the directory the function ends here:
	if(length(filelist) == 0){
		warning(paste("No pings files in the given directory ", event$event, sep=""))
		return()
	}
	# File numbers:
	if(length(filenr)==0 || identical(filenr, "all")){
		filenr <- seq_along(fileBaseNamesSansExt)
	}
	else if(!is.numeric(filenr)){
		warning("filenr must be numeric or \"all\"")
	}
	# Do not run more cores than pings, while keeping the possibly integer type of 'cores':
	if(length(filenr)==0){
		warning(paste("No pings files in the event", event$event))
	}
	else if(length(filenr)<cores){
		warning("The number of cores reduced to the number of pings-files")
		cores <- length(filenr)
	}

	# Order 'filenr' chronologically:
	indtStart <- lapply(indt[filenr], head, 1)
	lindtStart1 <- sapply(indtStart, length)==1
	if(!all(lindtStart1)){
		warning("Some files did not have any time steps, and were excluded from the compression")
		filenr <- filenr[lindtStart1]
		indtStart <- lapply(indt[filenr], head, 1)
	}
	filenr <- filenr[order(unlist(indtStart, use.names=FALSE))]

	# 'event_compr' is name of the directory for the TSD files (create it here and not before in case there are no pings files):
	if(!file.exists(event_compr)){
		suppressWarnings(dir.create(event_compr, recursive=TRUE))
	}
	
	# Create pings-filenames:
	pingsfiles <- getComprFileNames123(fileBaseNamesSansExt, event_compr, ext="pings")
	# Create beams-filenames:
	beamsfiles <- getComprFileNames123(fileBaseNamesSansExt, event_compr, ext="beams")
	# Create vessel-filenames:
	vesselfiles <- getComprFileNames123(fileBaseNamesSansExt, event_compr, ext="vessel")
	
	
	########### Processing: ###########
	##### Moving through the files of the event: #####
	# Declare the variable names:
	beamsnames <- labl.TSD("EKRaw2TSD_b")
	vesselnames <- labl.TSD("EKRaw2TSD_v")
	pingsnames <- labl.TSD("EKRaw2TSD_p")
	if(skipAngles){
		pingsnames <- setdiff(pingsnames, c("angl", "angt"))
	}
	
	## Check if all the files are present in the tsd-directory, in which case the function is terminated:
	#if(!ow){
	#	existingfiles <- list.files(event, recursive=TRUE, full.names=TRUE)
	#	if(length(existingfiles)){
	#		#cat("Existing files: \n",paste0(existingfiles, collapse="\n"), "\n")
	#		#cat("New files: \n",paste0(c(pingsfiles[[4]], beamsfile, vesselfile, ctdfile), collapse="\n"), "\n")
	#		ans <- readline("Overwrite existing files? (y/n)")
	#		if(!tolower(ans)=="y"){
	#			cat("Files already exist. Not overwriting\n")
	#			return()
	#		}
	#	}
	#}
	
	# Move through the list of pings files and read and possibly compress:
	indices <- papply(
		filenr, compr.event_oneFile_write, indt=indt, filelist=filelist, pingsfiles=pingsfiles, vesselfiles=vesselfiles, beamsfiles=beamsfiles, t="all", compress=compress, TIME=TIME, write=write, 
		# Inputs used in compr.TSD:
		tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, z0=z0, pingsnames=pingsnames, vesselnames=vesselnames, beamsnames=beamsnames, maxlenb=maxlenb, ..., cores=cores
	)
		
	
	
	# If there are files in the final merged directory, move the last file to the temporarily-merged directory (but only if it has 1 time step), and apply merging by time steps:
	x_final <- list.files(event_compr_final, full.names=TRUE)
	if(length(x_final)){
		# Split by file extension:
		ext_final <- file_ext(x_final)
		# Get the final pings, beams, and vessel files:
		x_final <- split(x_final, ext_final)
		names(x_final) <- sort(unique(ext_final))
		x_final <- x_final[names(x_final) %in% c("pings", "beams", "vessel")]
		# Check if the last files all have 1 time step, and move the files if so:
		lastfiles <- sapply(x_final, tail, 1)
		lastNumt <- read.TSDs(lastfiles, var="numt", clean=FALSE)
		lastNumt <- lastNumt[names(lastNumt)=="numt"]
		if(all(sapply(lastNumt, length)==1)){
			file.rename(lastfiles, file.path(dirname(x_final[1]), basename(lastfiles)))
		}
	}
	
	# Merge the last file of the current file number and the first file of the next:
	# Sort the files by name to keep the time ordered (changed on 2018-03-08):
	pingsfiles_filenr <- sort(unlist(lapply(pingsfiles[1:3], "[", filenr), use.names=FALSE))
	beamsfiles_filenr <- sort(unlist(lapply(beamsfiles[1:3], "[", filenr), use.names=FALSE))
	vesselfiles_filenr <- sort(unlist(lapply(vesselfiles[1:3], "[", filenr), use.names=FALSE))
	#pingsfiles <- merge_TSD_by_Time(pingsfiles_filenr, fun=weighted.mean)
	#beamsfiles <- merge_TSD_by_Time(beamsfiles_filenr, fun=head1)
	#vesselfiles <- merge_TSD_by_Time(vesselfiles_filenr, fun=head1)
	
	# 2018-09-12: Changed to apply the new function merge_TSD_by_var(), which must be run first with ind=NULL, and the using the output ind as input to do the merging:
	if(length(tres)>0){
		var <- "utim"
	}
	else if(length(xres)>0){
		var <- "sadv"
	}
	
	ind <- merge_TSD_by_var(vesselfiles_filenr, var=var)
	
	pingsfiles <- merge_TSD_by_var(pingsfiles_filenr, ind=ind, fun=weighted.mean)
	beamsfiles <- merge_TSD_by_var(beamsfiles_filenr, ind=ind, fun=head1)
	vesselfiles <- merge_TSD_by_var(vesselfiles_filenr, ind=ind, fun=head1)
	#indt <- seq_len(sum(read.TSDs(pingsfiles, var="numt", merge=TRUE, cores=cores)$numt))
	
	# Then merge the files using the linked option:
	#combine.TSD(pingsfiles, dir=event_compr_final, linked=list(beamsfiles, vesselfiles), adds=list(indt=indt), clear_along=TRUE)
	# combine.TSD(pingsfiles, dir=event_compr_final, linked=list(beamsfiles, vesselfiles), clear_along=clear_along, skipLast=TRUE) # WHY SKIPLAST?????????????????
	#combine.TSD(pingsfiles, dir=event_compr_final, linked=list(beamsfiles, vesselfiles), clear_along=clear_along)
	
	# 2018-09-19:
	#combine.TSD(vesselfiles, dir=event_compr_final, linked=list(beamsfiles, pingsfiles), clear_along=clear_along)
	combine.TSD(pingsfiles, dir=event_compr_final, linked=list(beamsfiles, vesselfiles), clear_along=clear_along)
	
	# Rename all files by removing "_1" just before the file extension:
	l <- list.files(event_compr_final, full.names=TRUE)
	file_path_sans_ext_l <- file_path_sans_ext(l)
	newl <- paste(substr(file_path_sans_ext_l, 1, nchar(file_path_sans_ext_l)-2), file_ext(l), sep=".")
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
#' (Internal) Read and write one Simrad raw file.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD write.TSD numt.TSD
#'
#' @keywords internal
#' @export
#' @rdname compr.event
#' 
compr.event_oneFile_write <- function(i, indt, filelist, pingsfiles, vesselfiles, beamsfiles, t="all", drop=TRUE, compress=FALSE, TIME=NULL, write=TRUE, 
	# Inputs used in compr.TSD:
	tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("mean","median"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, z0=0, pingsnames=NULL, vesselnames=NULL, beamsnames=NULL, dumpfiles=NULL, maxlenb=NULL, ...){

	# Function used for writing the first, last and the rest of the time steps in three separate files:
	write.TSD123 <- function(data, names, files, i, numt, ...){
		# Write the first ping:
		if(length(dumpfiles)){
			write(paste("TSD1_", i), dumpfiles[i], append=TRUE)
		}
		this = subset_TSD(data[names], ind=list(1), pad="start", drop=FALSE)
		TSD::write.TSD(this, files$comprFileNames1[i], numt=1, ...)
		# Write the middle pings:
		if(numt>1){
			if(length(dumpfiles)){
				write(paste("TSD2_", i), dumpfiles[i], append=TRUE)
			}
			thisindt <- seq(2, max(2,numt-1))
			this = subset_TSD(data[names], ind=list(thisindt), pad="start", drop=FALSE)
			TSD::write.TSD(this, files$comprFileNames2[i], numt=length(thisindt), ...)
		}
		# Write the last ping:
		if(numt>2){
			if(length(dumpfiles)){
				write(paste("TSD3_", i), dumpfiles[i], append=TRUE)
			}
			this = subset_TSD(data[names], ind=list(numt), pad="start", drop=FALSE)
			TSD::write.TSD(this, files$comprFileNames3[i], numt=1, ...)
		}
	}
		
	# Read the data from the raw file:
	#cat("Reading pings file number ", i, " of ", length(pingsfiles$filenames), " (", basename(pingsfiles$filenames[i]), ") ...\n", sep="")
	if(length(dumpfiles)){
		write(paste("File", i), dumpfiles[i], append=TRUE)
	}
		
	#data = read.event(filelist, t=indt[[i]], var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE, onestep=FALSE)
	data = read.event(filelist, t=indt[[i]], var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE, onestep=2)
	
	# Discard pings with too long beams:
	maxLenbPerPing = if(length(dim(data$lenb))==2) apply(data$lenb, 2, max) else max(data$lenb)
	if(length(maxlenb)>0 && any(maxLenbPerPing>maxlenb)){
		validPings = maxLenbPerPing <= maxlenb
		data = read.event(filelist, t=indt[[i]][validPings], var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE)
		write(paste("Discarded pings:", paste(which(!validPings), collapse=", ")), dumpfiles[i], append=TRUE)
	}
	#data = read.event(filelist, t=head(indt[[i]], 3), var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE)
	
	# Abort empry data:
	numt = numt.TSD(data)
	if(numt==0){
		warning(paste("No valid pings in file", i))
		return()
	}
	# Compress the data:
	if(compress){
		if(length(dumpfiles)){
			write(paste("Compr", i), dumpfiles[i], append=TRUE)
		}
		data <- compr.TSD(data, tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, z0=z0, ...)
		# Update numt for the compressed data:
		numt = numt.TSD(data)
	}
	
	if(length(dumpfiles)){
		write(paste("Write", i), dumpfiles[i], append=TRUE)
	}
	if(write){
		# Write the first time step, the last time step and the center time steps, so that the last and first time steps can be merged between data coming from consecutive raw files:
		if(compress){
			write.TSD123(data, pingsnames, pingsfiles, i, numt, header=list(dtyp=list(vbsc="floa")))
			write.TSD123(data, beamsnames, beamsfiles, i, numt)
			write.TSD123(data, vesselnames, vesselfiles, i, numt)
		}
		else{
			TSD::write.TSD(data[pingsnames], pingsfiles$filenames[i], numt=numt, header=list(dtyp=list(vbsc="floa")))
			TSD::write.TSD(data[beamsnames], beamsfiles$filenames[i], numt=numt)
			TSD::write.TSD(data[vesselnames], vesselfiles$filenames[i], numt=numt)
		}
	}
	if(length(dumpfiles)){
		write(paste("End", i), dumpfiles[i], append=TRUE)
	}
	#invisible(data)
	#rm(data)
	#gc()
	return(filelist[i])
}	
#*********************************************
#*********************************************
#' This funciton compresses the data stored in 'x' by averaging (either mean, trimmed mean or median) in depth and time bins. This is designed primarily for echosounder beams.
#'
#' @importFrom ccaPP fastMedian
#' @importFrom data.table := data.table key setkeyv
#' @importFrom SimradRaw soundbeam_range
#' @importFrom TSD dim_all global2car labl.TSD NAs strff utim.TSD utim2mtim zeros
#' @importFrom utils tail head
#'
#' @export
#' @rdname compr.event
#'
compr.TSD <- function(data=NULL, tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("mean","median"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, z0=0, drop=FALSE, ...){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-09-25 - Clean version.

	# Function stolen from http://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value. Repeats the last non-missing value. Used when filling in data in empty time intervals:
	repeatNonNA <- function(x){
		olddim = dim(x)
		ind = which(!is.na(x))
		# If it begins with a missing, add the first position to the indices:
		if(is.na(x[1])){
			ind = c(1, ind)
		}
		# Repeat the values at these indices diffing the indices + length yields how often :
		x = rep(x[ind], times=diff(c(ind, length(x) + 1) ))
		dim(x) = olddim
		x
	}		
	# Function used for applying the function 'funt':
	applyFunt <- function(x, tindex, funt, numt=NULL){
		applyFuntOne = function(xx){
			# Expand the indices to the array:
			olddim = dim_all(xx)
			len = length(xx)
			n1 = len/length(tindex)
			ind1 = rep(tindex, each=n1)
			ind2 = rep(seq_len(n1), length(tindex))
			
			# Apply the function 'funt' to the time intervals:
			out = c(tapply(xx, list(ind2, ind1), function(xxx) if(is.character(xxx)) head(xxx,1) else funt(xxx)))
			utindex = unique(tindex)
			dim(out) = c(olddim[-length(olddim)], length(utindex))
			
			# If numt is given larger than the number of unique time steps, return an expanded array with NAs at missing time intervals:
			if(length(numt) && numt>length(utindex)){
				outWithNA = NAs(olddim[-length(olddim)], numt)
				# Special care for one-dimensional vectors:
				if(length(olddim)==1){
					outWithNA[utindex] = out
				}
				else{
					outWithNA[, utindex] = out
				}
				out = repeatNonNA(outWithNA)
			}
			out
		}
		lapply(x, applyFuntOne)
	}
	
	if(length(data$vbsc)==0){
		return()
	}
	if(any(sapply(c(tres, xres, zres, rres, bres), length)>0)){
		# Merge the 'adds' and the data:
		if(length(adds)>0){
			data[names(adds)] = adds
		}
		data$utim = utim.TSD(data[labl.TSD("t")])

		# Get the dimensions of the vbsc
		dimvbsc = c(max(data$lenb), data$numb[1], length(data$utim))
		outdim = dimvbsc
		rindex = NULL
		tindex = NULL
		bindex = NULL
		origin = c(data$lonv[origin], data$latv[origin])
		if(is.null(data$psxv) || is.null(data$psyv)){
			data$psyv = global2car(cbind(c(data$lonv), c(data$latv)), origin=origin)
			data$psxv = data$psyv[, 1]
			data$psyv = data$psyv[, 2]
		}
		
		# Get the names of each category:
		beamsnames = setdiff(intersect(names(data), labl.TSD("b")), labl.TSD("t"))
		vesselnames = intersect(names(data), labl.TSD("v"))
		timenames = intersect(names(data), labl.TSD("t"))
		
		if(is.character(funvbsc)){
			if(strff("mean", funvbsc[1])){
				funvbsc = function(xx) sum(xx, na.rm=TRUE)/length(xx)
			}
			else if(strff("median", funvbsc[1])){
				funvbsc = fastMedian
			}
			else{
				funvbsc = get(funvbsc[1])
			}
		}
		
		if(is.character(funt)){
			if(strff("mean", funt[1])){
				funt = function(xx) sum(xx, na.rm=TRUE)/length(xx)
			}
			else if(strff("median", funt[1])){
				funt = fastMedian
			}
			else{
				funt = get(funt[1])
			}
		}
		
		
		###########################################
		########## 1. Compress in range: ##########
		###########################################
		# Compression in range 
		if(length(rres)>0){
			rindex = rep(soundbeam_range(data[c("lenb", "rres", "asps", "sint")], "mid"), prod(tail(dimvbsc,2)))
			# Split the data into range bins:
			if(length(rres)==1){
				if(length(data$rres)==0){
					# Get range resolution:
					data$rres = data$asps[1] * data$sint[1]/2
				}
				# Change rres and lenb to generate range intervals:
				rindexintervals = soundbeam_range(data[c("lenb", "rres", "asps", "sint")], pos="grid", adds=list(lenb=rep(ceiling(data$lenb[1] * data$rres[1] / rres), data$numb[1]), rres=rep(rres, l=data$numb[1])))
				#rangeintervals = seq(floor(min(rindex)/rres), ceiling(max(rindex)/rres)) * rres
			}
			else{
				rindexintervals = rres
			}
			rindex = findInterval(rindex, rindexintervals, all.inside=TRUE)
		}
		### 2. ... or compress in depth: ###
		else if(length(zres)>0){
			rindex = psx.TSD(data, pad=TRUE, split=split, ...)$psz
			# Split the data into depth bins:
			if(length(zres)==1){
				rindexintervals = seq(floor(min(rindex, na.rm=TRUE)/zres), z0) * zres
			}
			else{
				rindexintervals = zres
			}
			rindex = findInterval(rindex, rindexintervals, all.inside=TRUE)
			# Reverse the intervals, since depth is negative:
			rindex = length(rindexintervals) - rindex
			# Modify the transducer depth to z0:
			if(length(dim(data$psze))==2){
				data$psze = matrix(z0, nrow=1, ncol=dimvbsc[3])
			}
			else{
				data$psze = rep(z0, length.out=dimvbsc[3])
			}
			
		}
		# Alter the dimensions, and change the length of beams and range resolution:
		if(length(rindex)){
			outdim[1] = length(rindexintervals)-1
			# Change the length of the beams and the radial resolution
			data$lenb = matrix(length(rindexintervals)-1, nrow=dimvbsc[2], ncol=dimvbsc[3])
			data$rres = rep(tail(diff(rindexintervals),1), l=dimvbsc[3])
		}
		else{
			rindex = rep(seq_len(dimvbsc[1]), prod(dimvbsc[-1]))
		}
		
		##########################################
		########## 3. Compress in time: ##########
		##########################################
		if(length(tres)>0){
			vec = data$utim
			# Split the data into time bins:
			if(length(tres)==1){
				tindexintervals = seq(floor(min(vec, na.rm=TRUE)/tres), ceiling(max(vec, na.rm=TRUE)/tres)) * tres
			}
			tindex = findInterval(vec, tindexintervals, all.inside=TRUE)
		}
		### 4. ... or compress over sailed distance: ###
		else if(length(xres)>0){
			# Convert to meters:
			#tindex = data$utim[1] + (data$sadv - min(data$sadv)) * 1852
			# The above added time to distance, which does not make sense:
			vec = data$sadv * 1852
			# Split the data into sailed distance bins:
			if(length(xres)==1){
				tindexintervals = seq(floor(min(vec, na.rm=TRUE)/xres), ceiling(max(vec, na.rm=TRUE)/xres)) * xres
			}
			tindex = findInterval(vec, tindexintervals, all.inside=TRUE)
		}
		# Alter the dimensions, and subset the beams and vessel data:
		if(length(tindex)){
			
			# Add log, lon, lat and time start and end:
			start_t <- which(!duplicated(tindex))
			end_t <- c(start_t[-1], length(tindex))
			
			#data$start_t <- start_t
			#data$end_t <- end_t
			
			data$log1 <- data$sadv[start_t]
			data$lon1 <- data$lonv[start_t]
			data$lat1 <- data$latv[start_t]
			data$utm1 <- data$utim[start_t]
			
			data$log2 <- data$sadv[end_t]
			data$lon2 <- data$lonv[end_t]
			data$lat2 <- data$latv[end_t]
			data$utm2 <- data$utim[end_t]
			
			
			# Get the unique time steps indices:
			outdim[3] = length(tindexintervals) - 1
			
			# Extract beams data. Here the use of numt either fills the empty time intervals with NAs:
			if(length(dim(data$freq))==2){
				beamsToBeChanged = sapply(data[beamsnames], function(xx) if(length(xx)>0) tail(dim_all(xx),1)==dimvbsc[3] else FALSE)
				data[beamsnames][beamsToBeChanged] = applyFunt(data[beamsnames][beamsToBeChanged], tindex, funt, numt=outdim[3])
			}
			# Extract vessel data:
			thesevesselnames = setdiff(vesselnames, timenames)
			data[thesevesselnames] = applyFunt(data[thesevesselnames], tindex, funt, numt=outdim[3])
			#data[thesevesselnames] = lapply(data[thesevesselnames], function(xx) tapply(xx, tindex, funt))
			
							# Extract the time information using all intervals:
							diffTindexintervals = diff(tindexintervals)
							temp <- tindexintervals[-length(tindexintervals)] + diffTindexintervals/2
							if(length(tres)>0){
								data$utim <- temp
								#data$sadv <- 
							}
							### 4. ... or compress over sailed distance: ###
							else if(length(xres)>0){
								data$utim = (data$utm1 + data$utm1) / 2
								data$sadv <- temp
							}
							
							# Get the number of time steps in each compressed time bin:
							data$nmtc = zeros(outdim[3])
							data$nmtc[unique(tindex)] = table(tindex)
			
			#data$utim = tindexintervals[presentTindexintervals] + diffTindexintervals[presentTindexintervals]/2
			##data$utim = tindexintervals[-length(tindexintervals)] + diff(tindexintervals)/2
			
			# Add Matalb time (but remove this first and generate from the utim):
			data$mtim <- NULL
			data$mtim = mtim.TSD(data)
			data$indt = NULL
			#data$tnxi = tindexintervals
			#data$tndx = tindex
			tindex = rep(tindex, each=prod(head(dimvbsc,2)))
		}
		else{
			tindex = rep(seq_len(dimvbsc[3]), each=prod(dimvbsc[1:2]))
		}
		
		#############################################
		########## 5. Compress over beams: ##########
		#############################################
		oldindi = seq_len(data$numb[1])
		# We need this to separate between beams:
		bindex = rep(rep(oldindi, each=dimvbsc[1]), dimvbsc[3])
		if(length(bres)>0){
			# Split the data into beam bins:
			if(length(bres)==1){
				bindexintervals = 0.5 + seq(0, ceiling(max(oldindi)/bres)) * bres
			}
			else{
				bindexintervals = bres
			}
			# Set the new dira and dire:
			olddiredirainnew = findInterval(oldindi, bindexintervals, all.inside=TRUE)
			
			# Treat beams data given independent of time (equal for all pings):
			if(length(dim(data$dira))==0){
				beamsToBeChanged = sapply(data[beamsnames], function(xx) length(xx)==max(sapply(data[beamsnames], length)))
				data[beamsnames][beamsToBeChanged] = lapply(data[beamsnames][beamsToBeChanged], function(xx) by(xx, olddiredirainnew, funt))
				data$numb = length(data$dira)
				bindex = findInterval(bindex, bindexintervals, all.inside=TRUE)
				outdim[2] = length(bindexintervals)-1
			}
			# Treat beams data given per ping:
			else{
				beamsToBeChanged = sapply(data[beamsnames], function(xx) NROW(xx)==max(sapply(data[beamsnames], NROW)))
				data[beamsnames][beamsToBeChanged] = lapply(data[beamsnames][beamsToBeChanged], function(xx) do.call(rbind, by(xx, olddiredirainnew, function(yy) apply(yy, 2, funt))))
				data$numb = apply(data$dira, 2, length)
				dim(data$numb) <- c(1, length(data$numb))
				bindex = findInterval(bindex, bindexintervals, all.inside=TRUE)
				outdim[2] = length(bindexintervals)-1
			}
		}
		
		
		# Apply the compression
		if(!skipAngles && length(data$angt)){
			DT = data.table(vbsc=c(data$vbsc), angt=c(data$angt), angl=c(data$angl))
		}
		else{
			data$angt = NULL
			data$angl = NULL
			DT = data.table(vbsc=c(data$vbsc))
		}
		
		# Add the indices:
		presentKeys = NULL
		if(length(tindex)){
			DT[,tindex:=tindex]
			presentKeys = c(presentKeys, "tindex")
		}
		if(length(bindex)){
			DT[,bindex:=bindex]
			presentKeys = c(presentKeys, "bindex")
		}
		if(length(rindex)){
			DT[,rindex:=rindex]
			presentKeys = c(presentKeys, "rindex")
		}
		
		
		setkeyv(DT, presentKeys)
		
		data$vbsc = NAs(outdim)
		temp = DT[, funvbsc(vbsc), by=key(DT)]
		temp2 = as.matrix(temp[, rev(presentKeys), with=FALSE])
		data$vbsc[temp2] = temp$V1
		
		
		if(!skipAngles && length(data$angt)){
			warning("Electric angles are deprecated")
			
			data$angt = NAs(outdim)
			temp = DT[, funvbsc(angt), by=key(DT)]
			data$angt[as.matrix(temp[, rev(presentKeys), with=FALSE])] = temp$V1
			
			data$angl = NAs(outdim)
			temp = DT[, funvbsc(angl), by=key(DT)]
			data$angl[as.matrix(temp[, rev(presentKeys), with=FALSE])] = temp$V1
		}
		
		data[c("psxx", "psyx", "pszx")] = psx.TSD(data, pad=TRUE, split=split, ...)
	}
	else{
		warning("No compression")
	}
	
	if(drop){
		data$vbsc <- drop(data$vbsc)
	}
		
	
	##### Output #####
	invisible(data)
}
