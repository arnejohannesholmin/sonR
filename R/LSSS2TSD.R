#*********************************************
#*********************************************
#' Convert from LSSS Profos files to TSD files.
#'
#' \code{LSSS2TSD} Converts from LSSS Profos export files to TSD files, readable by R using the functions read.TSD, read.TSDs and read.event.
#'
# See EKRaw2TSD for parameter descriptions.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom data.table rbindlist
#' @importFrom SimradRaw soundbeam_range
#' @importFrom TSD arr.ind2ind write.TSD NAs zeros utim2ftim combine.TSD papply
#' @importFrom utils tail head read.table
#'
#' @export
#' @rdname EKRaw2TSD
#'
LSSS2TSD <- function(event=NULL, filesize=3e8, chunksize=1e8, keep.temp=FALSE, cores=1, ...){
	
	############### LOG: ###############
	# Start: 2016-12-06 - Clean version.

	##################################################
	##################################################
	########## Preparation ##########
	### Functions: ###
	# Get the string between two other strings:
	getStringBetween <- function(x, s1, s2){
		gsub(paste0(".*", s1, "\\s*|", s2, ".*"), "", x)
	}
	# Get the esnm from a string:
	getESNM <- function(x){
		toupper(getStringBetween(x, s1="Ping_", s2="-D"))
	}
	# Read sample data using read.table:
	readSampleData <- function(x){
		# Get sample data:
		read.table(x, header=TRUE, sep=",")
	}
	# Get the radial resolutions of a school file:
	readRange <- function(x){
		# Get sample data:
		read.table(x, header=TRUE, sep=",")$Range
	}
	# Get school ID with file identifiers as names:
	getSchoolID <- function(l, n=10){
		getSchoolIDOne <- function(x, n){
			l <- readLines(x, n=n)
			hash <- startsWith(l, "#")
			# Get school ID:
			schoolIDString <- "# School id: "
			schoolIDLine <- l[startsWith(l, schoolIDString)]
			as.numeric(substring(schoolIDLine, nchar(schoolIDString)+1))
		}
		getFileIDOne <- function(x, n){
			l <- readLines(x, n=n)
			hash <- startsWith(l, "#")
			# Get school ID:
			fileIDString <- "# Ping: "
			fileIDLine <- l[startsWith(l, fileIDString)]
			substring(fileIDLine, nchar(fileIDString)+1)
		}
		schoolID <- unlist(lapply(l, getSchoolIDOne, n=n))
		fileID <- unlist(lapply(l, getFileIDOne, n=n))
		names(schoolID) <- fileID
		schoolID
	}
	# Converts full data to segmentation data:
	#getSeg <- function(xx){
	#	present <- which(!is.na(xx))
	#	cbind(present, xx[present])
	#}
	# Reads acoustic data, range indices and beam indices:
	getSeg <- function(file, beams){
		x <- read.table(file, header=TRUE, sep=",")
		#out <- NAs(max(beams$lenb), beams$numb[1])
		ranges <- soundbeam_range(beams)
		edges <- soundbeam_range(beams, pos="edge")
		rangeInd <- findInterval(x$Range, edges)
		beamsInd <- colnames(x)
		beamsInd <- beamsInd[startsWith(beamsInd, "Beam")]
		beamsInd <- as.numeric(substring(beamsInd, 5)) + 1
		#out[rangeInd, beamsInd] <- as.matrix(x[,-1])
		ind <- arr.ind2ind(as.matrix(expand.grid(rangeInd, beamsInd)), c(max(beams$lenb), beams$numb[1]))
		data <- x[,-1]
		noNA <- !is.na(data)
		list(ind=ind[noNA], data=data[noNA])
	}
	# Read all the files for one ping, merge the schools, and write the pingsfile:
	getOnePing <- function(pingInd, files, outfiles, vessel, beams){
		# Read the segmentation masks:
		temp <- lapply(files[[pingInd]], getSeg, beams=beams)
		# Combine all matrices of indices and data:
		temp <- data.table::rbindlist(temp)
		# Sum over all indices (allowing for several overlapping schools, which is unlikely, but still the method is simple):
		temp <- tapply(temp[[2]], temp[[1]], sum)
		
		# Save compressed data (indices and segmented data):
		pings <- list()
		pings$vbsc = 10^(temp/10)
		pings$vxIX = as.numeric(names(temp))
		
		## Merge the data, counting NAs as 0:
		#out <- NAs(max(beams$lenb), beams$numb[1])
		#out[as.numeric(names(temp))] <- temp
		
		pings$utim <- vessel$utim[pingInd]
		#pings$vbsc <- 10^(out/10)
		pings$lenb <- beams$lenb[,pingInd]
		pings$numb <- beams$numb[pingInd]
		
		# Write the data to individual pings files and merge later:
		write.TSD(pings, con=outfiles[pingInd], numt=1)
	}
	### End of functions: ###
	
	
	# Default beams values:
	defaultBeams <- list(numb=NULL, lenb=NULL, psze=-7, diraSpan=c(0, 360), asps=1480, bwtx=NA, bwty=NA, eqba=NA, freq=NA, sacr=NA, tpow=NA, gai1=NA, gai2=NA, bmmd=NA)
	if(length(event)==0){
		return(defaultBeams)
	}
	
	ll <- list(...)
	beamsNames <- intersect(names(defaultBeams), names(ll))
	if(length(beamsNames)){
		defaultBeams[beamsNames] <- ll[beamsNames]
	}
	# Extract the defaults:
	numb <- defaultBeams$numb
	lenb <- defaultBeams$lenb
	psze <- defaultBeams$psze
	diraSpan <- defaultBeams$diraSpan
	if(any(diraSpan > 2*pi)){
		diraSpan <- diraSpan * pi/180
	}
	asps <- defaultBeams$asps
	bwtx <- defaultBeams$bwtx
	bwty <- defaultBeams$bwty
	eqba <- defaultBeams$eqba
	freq <- defaultBeams$freq
	sacr <- defaultBeams$sacr
	tpow <- defaultBeams$tpow
	gai1 <- defaultBeams$gai1
	gai2 <- defaultBeams$gai2
	bmmd <- defaultBeams$bmmd
	
	# Get the list of files, slit into the sub directories:
	dirs = list.dirs(event)
	l = list.files(event, recursive=TRUE, full.names=TRUE)
	lclean <- gsub(paste0(event, "/"), "", l, fixed=TRUE)
	lclean <- split(lclean, dirname(lclean))
	# Add the event again:
	lclean <- lapply(lclean, function(xx) file.path(event, xx))
	if(any(sapply(lclean[c("Export/ProfosSchoolSamples", "Reports/profos")], length)==0)){
		warning("The event must contain the directories \"Export/ProfosSchoolSamples\" and \"Reports/profos\"")
	}
	
	# Get the type of sonar:
	eventname = basename(event)
	esnm <- getESNM(head(lclean[["Export/ProfosSchoolSamples"]], 1))
	
	# Create the TSD directory and define file names:
	tsdDir <- file.path(event, esnm, "tsd")
	tsdDirIndividualPings <- file.path(event, esnm, "tsd_individual")
	suppressWarnings(dir.create(tsdDir, recursive=TRUE))
	suppressWarnings(dir.create(tsdDirIndividualPings, recursive=TRUE))
	
	vesselfile <- file.path(tsdDir, paste0(eventname, ".vessel"))
	beamsfile <- file.path(tsdDir, paste0(eventname, ".beams"))
	
	
	########### Processing: ###########
	# Define list storing the data:
	beams <- list()
	vessel <- list()
	pings <- list()
	
	
	
	# Read the compressed school data containing vessel information. This is the files ending with pp.txt (denoting "per ping"):
	ppFile <- lclean[["Reports/profos"]]
	ppFile <- ppFile[substring(ppFile, nchar(ppFile)-5)=="pp.txt"]
	ppData <- read.table(ppFile, header=TRUE)
	# Get school IDs:
	schoolID <- ppData$Id
	# Get the school IDs in the sample data files:
	schoolIDFiles <- getSchoolID(lclean[["Export/ProfosSchoolSamples"]])
	# Order the school files by the order in the ppFile:
	lclean[["Export/ProfosSchoolSamples"]] <- lclean[["Export/ProfosSchoolSamples"]][match(as.character(ppData$File), names(schoolIDFiles))]
	
	
	##### Vessel data: #####
	# 1. Get the time:
	vessel$utim <- unclass(as.POSIXct(paste(ppData$Date, ppData$Time), tz="UTC"))
	# The vessel times can be duplicated if there are more than one school in one or more pings:
	
	
	numt <- length(vessel$utim)
	# 2. Get the vessel position:
	vessel$lonv <- ppData$Ship.lon
	vessel$latv <- ppData$Ship.lat
	# 3. Get the vessel speed:
	vessel$ispv <- ppData$Ship.speed
	# 4. Get the vessel pitch (rtxv), roll (rtyv), heading (rtzv):
	vessel$rtxv <- NAs(numt)
	vessel$rtyv <- NAs(numt)
	vessel$rtzv <- (-ppData$Ship.heading * pi / 180) %% (2*pi)
	# 5. Get the vessel heave:
	vessel$pszv <- NAs(numt)
	# 6. Get the time offset of the vessel information:
	vessel$terr <- zeros(numt)
	
	# Sort the vessel information by the time steps:
	vessel <- lapply(vessel, function(xx) xx[order(vessel$utim)])
	
	# Group the school files into time steps:
	lclean[["Export/ProfosSchoolSamples"]] <- split(lclean[["Export/ProfosSchoolSamples"]], sort(unique(vessel$utim)))
	pingInd <- seq_along(lclean[["Export/ProfosSchoolSamples"]])
	
	# Uniquify the vesse data:
	duputim <- duplicated(vessel$utim)
	vessel <- lapply(vessel, "[", !duputim)
	
	
	##### Beams data: #####
	# 1. Get the sonar name:
	beams$esnm <- rep(esnm, numt)
	# 2. Get the UNIX time:
	beams$utim <- vessel$utim
	beams$asps <- rep(asps, numt)
	# Get number of beams:
	if(length(lenb)==0){
		if(esnm=="SN90"){
			numb <- 256
		}
		else if(esnm %in% c("SU90", "SX90", "SH90", "SX80", "SH80")){
			numb <- 64
		}
	}
	beams$numb <- rep(numb, numt)
	if(length(freq)==1 && is.na(freq)){
		warning("freq must be set manually")
	}
	beams$freq <- matrix(freq, nrow=numb, ncol=numt)
	beams$absr <- NAs(numb, numt)
	# Read the ranges of all files, and get the radial resolution as exactly as possible:
	dd <- papply(lclean[["Export/ProfosSchoolSamples"]], readRange, cores=cores, outfile=FALSE)
	beams$rres <- rep(mean(unlist(lapply(dd, diff))), numt)
	maxRange <- max(unlist(lapply(dd, max)))
	
	if(length(lenb)==0){
		lenb <- ceiling(maxRange / beams$rres[1])
	}
	beams$lenb <- matrix(lenb, nrow=numb, ncol=numt)
	
	# Get also the data from the first file of school data:
	beams$sint <- rep(round(2 * beams$rres[1] / asps * 1e6) * 1e-6, numt)
	# Update asps:
	beams$asps <- rep(2 * beams$rres[1] / beams$sint[1], numt)
	beams$plsl <- NAs(numt)
	beams$psze <- NAs(numt)
	if(length(lenb)==1 && is.na(lenb)){
		warning("lenb must be set manually")
	}
	# Create a sequence from diraSpan (if given as 0 to 2*pi, assume this is from beam 1 to beam 1):
	full <- identical(range(diraSpan), c(0,2*pi))
	deltadira <- diff(diraSpan) / (beams$numb[1] + if(full) 1 else 0)
	beams$dira <- matrix(seq(diraSpan[1], by=deltadira, l=numb), nrow=numb, ncol=numt)
	# Add 90 degrees and convert to radians
	beams$dire <- matrix((-ppData$Trans.tilt + 90) * pi/180, nrow=numb, ncol=numt, byrow=TRUE)
	# 1. Get the time:
	if((length(bwtx)==1 && is.na(bwtx)) || (length(bwty)==1 && is.na(bwty))){
		warning("bwtx and bwty must be set manually")
	}
	beams$bwtx <- zeros(numb, numt) + bwtx
	beams$bwty <- zeros(numb, numt) + bwty
	beams$eqba <- zeros(numb, numt) + eqba
	beams$sacr <- zeros(numb, numt) + sacr
	beams$tpow <- zeros(numb, numt) + tpow
	beams$gai1 <- zeros(numb, numt) + gai1
	beams$gai2 <- zeros(numb, numt) + gai2
	beams$bmmd <- zeros(numt) + bmmd

	# Write the vessel data:
	write.TSD(vessel, con=vesselfile, numt=numt)
	# Write the vessel data:
	write.TSD(beams, con=beamsfile, numt=numt)
	
	
	##### Pings: #####
	# Get file names of the pings:
	pingsfilesDateTime <- utim2ftim(as.numeric(unlist(names(lclean[["Export/ProfosSchoolSamples"]]))), "Dyyyymmdd-THHMMSS")
	pingsfiles <- file.path(tsdDirIndividualPings, paste0(eventname, pingsfilesDateTime, ".pings"))
	# Parallel converting of the school files:
	papply(pingInd, getOnePing, files=lclean[["Export/ProfosSchoolSamples"]], outfiles=pingsfiles, vessel=vessel, beams=beams, cores=cores, outfile=FALSE)
	
	# Simply merge the pings files (pad=FALSE ensures that segmented data, using the vxIX, are saved in a list, and not padded with NAs):
	pingsfiles = combine.TSD(pingsfiles, dir=tsdDir, reserve=TRUE, filesize=filesize, chunksize=chunksize, cores=cores, pad=FALSE, ...)
	
	# Delete the individual pings:
	if(!keep.temp){
		unlink(tsdDirIndividualPings, recursive=TRUE)
	}
	
	# Return individual lists or merge to one list:
	invisible(list(vesselfile=vesselfile, beamsfile=beamsfile, pingsfiles=pingsfiles))
	##################################################
	##################################################
	}
