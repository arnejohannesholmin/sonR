#*********************************************
#*********************************************
#' Identifies whether the acoustic system is a sonar, fishery sonar or echosounder, for exactly one ping is 'data' is given.
#'
#' @param esnm  the four character identifyer of the acoustical system. Currently implemented: echosounders: "EK60", "ME70"; sonars: "MS70", "SX80", "SH80", "SX90", "SH90".
#' @param bydirs  is TRUE if the elevation angle of the beams should be used to identify sonars by the requirement that none of the beams should be within 'margin' of the vertical orientation.
#' @param data  is an optional list of beam configuration data.
#' @param margin  is the margin around the vertical orientation outside which all beams must be to be identified as sonar.
#' @param bmmd  is the desired mean mode of the system, which should be set for instance for fishery sonars.
#' @param fishery  is TRUE to return TRUE only for fishery sonars (not including MS70).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname is.sonar
#'
is.sonar<-function(esnm="MS70", bydirs=FALSE, data=list(), margin=10*pi/180, bmmd=NULL, fishery=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-08-28 - Clean version.
	########### DESCRIPTION: ###########
	# Identifies whether the acoustic system is a sonar, fishery sonar or echosounder, for exactly one ping is 'data' is given.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---esnm--- the four character identifyer of the acoustical system. Currently implemented: echosounders: "EK60", "ME70"; sonars: "MS70", "SX80", "SH80", "SX90", "SH90".
	# ---bydirs--- is TRUE if the elevation angle of the beams should be used to identify sonars by the requirement that none of the beams should be within 'margin' of the vertical orientation.
	# ---data--- is an optional list of beam configuration data.
	# ---margin--- is the margin around the vertical orientation outside which all beams must be to be identified as sonar.
	# ---bmmd--- is the desired mean mode of the system, which should be set for instance for fishery sonars.
	# ---fishery--- is TRUE to return TRUE only for fishery sonars (not including MS70).
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Define echosounders:
	echsounders = c("ek60", "me70")
	# Define sonars:
	fishsonars = c("sx80", "sh80", "su80", "sx90", "sh90", "su90")
	sonars = c("ms70", fishsonars)
	# Test:
	if(is.list(esnm)){
		data = esnm
		esnm = esnm$esnm[1]
		}
	issonar = if(fishery) tolower(esnm) %in% fishsonars else tolower(esnm) %in% sonars
	#isechsounder = tolower(esnm) %in% echsounders
	
	
	########## Execution and output ##########
	#if(bydirs==FALSE && any(issonar,isechsounder)){
	if(!bydirs){
		issonar
		}
	else if(length(data)>0){
		if(length(unique(data$bmmd))>1){
			if(length(bmmd)==0){
				stop("The beam mode 'bmmd' must be given in the beam configuration data are given for several beam modes (length(unique(data$bmmd))>1)")
				}
			all(data$dire[data$bmmd==bmmd]<pi-margin)
			}
		else{
			all(data$dire<pi-margin)
			}
		}
	##################################################
	##################################################
	}
