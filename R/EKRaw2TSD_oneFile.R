#*********************************************
#*********************************************
#' (Internal) Read and Simrad raw file.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw NMEA2vessel readEKRaw readEKRaw_power2sv.TSD readEKRaw_stripNA getRangeOffsetInUnitsOfSamples
#' @importFrom TSD listOfEqual2array mergeListKeepDimensions NAs
#'
#' @importFrom utils tail head
#' @importFrom SimradRaw offsetSamples
#'
#' @export
#' @rdname EKRaw2TSD_oneFile
#' 
EKRaw2TSD_oneFile <- function(i, filelist,  prenumt=10, t="all", endian="little", timeOffset=0, minTimeDiff=Inf, msg=TRUE, na.rm=TRUE, correctTime=TRUE, bmmd=NULL, TVG.exp=2, dira_offset=0, cali=TRUE, toTS=FALSE, psze=-7, skipAngles=TRUE, cleanNMEA=1, apply.range.offset=FALSE, thr1m=FALSE){

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
	
	# If specified, split the channels into separate time steps with beam modes. This is only applied if 'bmmd' is given ???????? When is this ever used ??????? (2016-12-12):
	if(length(bmmd) == thisd$header$transceivercount){
		head2 <- function(x, n=1, N=2){
			if(NCOL(x)==N){
				head(t(x), n)
			}
			else{
				head(x, n)
			}
		}
		numb <- thisd$header$transceivercount
		thisd <- lapply(thisd$data, function(x1) lapply(x1, function(x2) if(is.list(x2)) lapply(x2, head2, 1, numb) else head2(x2, 1, numb)))
	}
	
	########## (2) Store beams variables: ##########
	numt <- NCOL(thisd$data$pings$time)
	# (2.1) esnm - the name of the system:
	beams$esnm <- EKRaw2TSD_getesnm(thisd, numt)
	# (2.3) mtim - Matlab time:
	beams$mtim <- thisd$data$pings$time[1,]
	# (2.4) asps - average speed of sound:
	beams$asps <- thisd$data$pings$soundvelocity[1,]
	# (2.5) numb - number of beams:
	beams$numb <- rep(thisd$header$transceivercount,numt)
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
	if(length(thisd$data$pings$transducerdepth) && all(thisd$data$pings$transducerdepth!=0)){
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
		valid <- !is.na(pings$mtim)
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
