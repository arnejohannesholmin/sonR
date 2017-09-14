#*********************************************
#*********************************************
#' Reads the time and position of a CTD file in the cnv format.
#'
#' @param x  is a CTD file in the cnv format.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname CTDtimeandpos
#'
CTDtimeandpos<-function(x){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-10-31 - Clean version.
	########### DESCRIPTION: ###########
	# Reads the time and position of a CTD file in the cnv format.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is a CTD file in the cnv format.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	if(isTRUE(file.exists(x))){
		x=readLines(x)
		}
	start = substr(x,1,18)
	atlat = start == "* NMEA Latitude = "
	atlon  = start == "* NMEA Longitude ="
	attime  = start == "* System UpLoad Ti"
	utim = unclass(as.POSIXct(strptime(substring(x[attime],24),"%b %d %Y %H:%M:%S")))
	lat = strsplit(substr(x[atlat],19,nchar(x[atlat])-1)," ",fixed=TRUE)[[1]]
	lat = as.numeric(lat[nchar(lat)>0])
	lon = strsplit(substr(x[atlon],19,nchar(x[atlon])-1)," ",fixed=TRUE)[[1]]
	lon = as.numeric(lon[nchar(lon)>0])
	lat = lat[1] + lat[2]/60
	lon = lon[1] + lon[2]/60
	North = substring(x[atlat],nchar(x[atlat])) == "N"
	East = substring(x[atlon],nchar(x[atlon])) == "E"
	if(!North){
		lat = -lat
		}
	if(!East){
		lon = -lon
		}
	
	########## Execution and output ##########
	c(utim=utim, lat=lat, lon=lon)
	##################################################
	##################################################
	}
