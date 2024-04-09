#*********************************************
#*********************************************
#' Converts a CTD file to the TSD format for the variables relevant to echoIBM().
#'
#' @param f  is the CTD file in cnv format.
#' @param outfile  is file to which the CTD data are written.
#' @param hpr0  is the pressure at sea level, which must be given externally. Defaulted to 101325 pascal.
#' @param Pain  is TRUE if pressure "ihpr" is given in Pascal and FALSE if given in decibar relative to surface pressure (10000 Pascal, giving values approximately equivalent to water depth).
#' @param vessel  A list of vessel info by which the closest CTD station is found as the minimum of the sum of distance in meters and seconds.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD global2car mtim2utim write.TSD
#' @importFrom tools file_ext
#' @importFrom utils tail head
#'
#' @export
#' @rdname CTD2TSD
#'
CTD2TSD <- function(f=NULL, outfile=NULL, hpr0=NULL, Pain=FALSE, vessel=list()){
	
	############### LOG: ###############
	# Start: 2013-11-16 - Clean version.
	
	##### Preparation #####
	if(length(f) == 0){
		f <- file.path(Acoustics_datasets_directory(), "CTDs")
	}
	if(length(f) == 1){
		if(isTRUE(file.info(f)$isdir)){
			f <- list.files(f, recursive=TRUE, full.names=TRUE)
		}
	}
	f <- f[file_ext(f) == "cnv"]
	if(length(f) == 0){
		return()
	}
	
	# Get the CTD file times and positions:
	CTD_info <- matrix(unlist(lapply(f, CTDtimeandpos)), byrow=TRUE, ncol=3)
	dist <- sqrt(rowSums(global2car(CTD_info[,2:3], origin=c(vessel$lonv[1], vessel$latv[1], 0))^2))
	closest <- which.min(dist + abs(CTD_info[,1] - utim.TSD(vessel)[1]))
	f <- f[closest]
	
	if(isTRUE(file.info(f)$isdir)){
		f <- list.files(f, full.names=TRUE)
		ext <- file_ext(f)
		f <- f[ext == "cnv"]
		if(ext != "cnv"){
			warning("The input directory must contain a cnv file")
			return(list())
		}
	}
	# Read the data file:
	ext <- file_ext(f)
	if(ext!="cnv"){
		warning("The input file must be a cnv file")
		return(list())
	}
	l <- readLines(f)
	
	# Define strings which recognize the variables:
	namesstr <- "# name"
	endstr <- "*END*"
	atnames <- grep(namesstr, l, fixed=TRUE)
	atend <- grep(endstr, l, fixed=TRUE)
	
	
	##### Execution and output #####
	# Get the UNIX time, latitude, longitude:
	utimlatlon <- CTDtimeandpos(l)
	utim <- utimlatlon[1]
	latc <- utimlatlon[2]
	lonc <- utimlatlon[3]
	
	# Get the variable names
	atequal <- sapply(gregexpr("=", l[atnames], fixed=TRUE), head, 1)
	atcolon <- sapply(gregexpr(":", l[atnames], fixed=TRUE), head, 1)
	varnames <- substr(l[atnames], atequal+2, atcolon-1)
	
	# Read the data:
	atdata <- seq(atend+1, length(l)-1)
	data <- strsplit(l[atdata]," ")
	data <- lapply(data, function(x) x[nchar(x)>0])
	data <- unlist(data)
	dim(data) <- c(length(varnames), length(atdata))
	data <- t(data)
	colnames(data) <- varnames
	
	### Merge the two sensors: ###
	# Temperature:
	temppresent <- which(varnames %in% c("t068C","t168C"))
	temppresent <- matrix(as.numeric(data[,temppresent,drop=FALSE]), ncol=length(temppresent))
	if(length(temppresent)>0){
		temp <- rowMeans(temppresent)
	}
	else{
		temp <- NULL
	}
	
	# Salinity:
	sltypresent <- which(varnames %in% c("sal00","sal11"))
	sltypresent <- matrix(as.numeric(data[,sltypresent,drop=FALSE]), ncol=length(sltypresent))
	if(length(sltypresent)>0){
		slty <- rowMeans(sltypresent)
	}
	else{
		slty <- NULL
	}
	
	# Pressure:
	ihprpresent <- which(varnames %in% c("prdM", "prDM"))
	ihprpresent <- matrix(as.numeric(data[,ihprpresent,drop=FALSE]), ncol=length(ihprpresent))
	if(length(ihprpresent)>0){
		ihpr <- rowMeans(ihprpresent)
	}
	else{
		ihpr <- NULL
	}
	
	# Sound speed:
	ispspresent <- which(varnames %in% c("svCM"))
	ispspresent <- matrix(as.numeric(data[,ispspresent,drop=FALSE]), ncol=length(ispspresent))
	if(length(ispspresent)>0){
		isps <- rowMeans(ispspresent)
	}
	else{
		isps <- NULL
	}
	
	# Gravitational constant:
	gacc <- 9.780327 * (1 + 0.0053024*sin(latc)^2 - 0.0000058*sin(2*latc)^2)
	
	# Sound speed:
	rho0present <- which(varnames %in% c("sigma-t00"))
	rho0present <- matrix(as.numeric(data[,rho0present, drop=FALSE]),ncol=length(rho0present))
	if(length(rho0present)>0){
		rho0 <- rowMeans(rho0present)+1000
	}
	else{
		rho0 <- NULL
	}
		
	# Air pressure default (information exists in reflog files, but these are not found at the moment):
	if(length(hpr0)==0){
		hpr0 <- 101325/10000
	}
	
	# Get the depths:
	pszc <- getzfromctd(list(ihpr=ihpr, rho0=rho0, gacc=gacc, hpr0=hpr0),Pain=Pain)
	
	# The output:
	out <- list(utim=utim, lonc=lonc, latc=latc, pszc=pszc, ihpr=ihpr, temp=temp, slty=slty, isps=isps, rho0=rho0, gacc=gacc, hpr0=hpr0)
		
	# Write the TSD file:
	if(length(outfile)>0){
		write.TSD(out, outfile,numt=1)
	}
	
	# Return:
	out
}
