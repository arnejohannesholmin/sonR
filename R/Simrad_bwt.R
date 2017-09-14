#*********************************************
#*********************************************
#' Converts from Simrad beam widths to x and y beam widths.
#'
#' @param data  is a list containing the required variables 'dirx' (vector of x-beam-directions), 'diry' (vector of y-beam-directions), 'dirz' (vector of z-beam-directions), and 'esnm' (name of the echosounder/sonar).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname Simrad_bwt
#'
Simrad_bwt<-function(data){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-10-31 - Clean version.
	########### DESCRIPTION: ###########
	# Converts from Simrad beam widths to x and y beam widths.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the required variables 'dirx' (vector of x-beam-directions), 'diry' (vector of y-beam-directions), 'dirz' (vector of z-beam-directions), and 'esnm' (name of the echosounder/sonar).
	
	
	##################################################
	##################################################
	### 'bwtx' is the horizontal beam width and 'bwty' is the vertical beam width. For vertical beams, 'bwtx' is the beam width atwarth the ship and 'bwty' is the beam width along the ship.
	# The interpretation of the terms "alongship" and "atwarthship" differ for EK60 and ME70 (intuitive interpretation) and MS70 (counterintuitive interpretation):
	# For the EK60, the x-axis points atwarthship (straigth out on the starboard side of the vessel), and alongship is along the ship:	
	if(sonR_implemented(data, "SBE")){
		if(length(data$bwtx)==0){
			data$bwtx=data$bwtt * pi/180
			}
		if(length(data$bwty)==0){
			data$bwty=data$bwtl * pi/180
			}
		}
	# ?????For the ME70, the x-axis points alongship. For the beams of the ME70 that are on the starboard side of the vessel, the coordinate system of the vessel is first rotated by -90 degrees around the z axis, and then by an angle > 90 around the x axis, so that the x axis points along the ship backwards:?????
	else if(sonR_implemented(data, "MBE")){	
		if(length(data$bwtx)==0){
			data$bwtx=data$bwtl * pi/180
			}
		if(length(data$bwty)==0){
			data$bwty=data$bwtt * pi/180
			}
		}
	# For the MS70, the "0"-mark is not in the direction of the vessel, but pointing upwards on the port side of the vessel. Thus the atwarthship angle is actually along ship, corresponding to 'bwtx':
	else if(sonR_implemented(data, "MBS")){	
		if(length(data$bwtx)==0){
			data$bwtx=data$bwtt * pi/180
			}
		if(length(data$bwty)==0){
			data$bwty=data$bwtl * pi/180
			}
		}
	else if(sonR_implemented(data, "OFS")){	
		#stop("sh80, sx80, sh90 or sx90 not yet implemented")
		# This needs to be checked:
		if(length(data$bwtx)==0){
			data$bwtx=data$bwtl * pi/180
			}
		if(length(data$bwty)==0){
			data$bwty=data$bwtt * pi/180
			}
		}
	else{
		warning("Unknown sonar or echosounder. 'bwtx' chosen to be 'bwtl', and 'bwty' chosen to be 'bwtt'")
		if(length(data$bwtx)==0){
			data$bwtx=data$bwtl * pi/180
			}
		if(length(data$bwty)==0){
			data$bwty=data$bwtt * pi/180
			}
		}
	data
	##################################################
	##################################################
	}
