#*********************************************
#*********************************************
#' Converts from Simrad beam directions to azimuth and elevation angles.
#'
#' @param data  is a list containing the rquired variables 'dirx' (vector of x-beam-directions), 'diry' (vector of y-beam-directions), 'dirz' (vector of z-beam-directions), and 'esnm' (name of the echosounder/sonar).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname Simrad_dir
#'
Simrad_dir<-function(data, dira_offset=0){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-10-31 - Clean version.
	# Last: 2015-03-17 - Added support for the old variables 'dirl' and 'dirt', and changed name from "Simrad_dirl_dirt2dira_dire" to "Simrad_dir".
	########### DESCRIPTION: ###########
	# Converts from Simrad beam directions to azimuth and elevation angles.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the rquired variables 'dirx' (vector of x-beam-directions), 'diry' (vector of y-beam-directions), 'dirz' (vector of z-beam-directions), and 'esnm' (name of the echosounder/sonar).
	
	
	##################################################
	##################################################
	# If the old variables 'dirl' and 'dirt' are the only direction angles in 'data', obtain 'dirx' and 'diry' from these:
	if(sum(length(data$dirx)>0, length(data$diry)>0)<2 && sum(length(data$dirl)>0, length(data$dirt)>0)==2){
		# <<< COPIED FROM THE OLD FUNCTION read.event_get.dira_dire >>>
		# The interpretation of "along" and "athwart" ship differs for echosounders and sonars, where the sonars that we know of are defined with angles representing the elevation of the angles given in the 'dirx' field. For a sonar looking typically sideways, this would be athwart the ship. For echosounders the interpretations are reversed:
		if(sonR_implemented(data, c("MBS","OFS"))){
			data$diry = data$dirl
			data$dirx = data$dirt
			}
		else{
			data$dirx = data$dirl
			data$diry = data$dirt
			}
		}
	if(sum(length(data$dira)>0, length(data$dire)>0)<2){
		# For the MS70 sonar, the direction angles are in spherical coordinates, but given so that data$dirx ("along ship", the radian version of diry) is the azimuth angle deviation from the negative x-axis, and data$diry ("atwarth ship" the radian version of dirx) is the elevation angle deviation from the x-axis. Here the terms "along ship" and "atwarth ship" are used in the intuitive sense, and not according to the definition of the MS70 sonar. 
		if(sonR_implemented(data, "MBS")){
			data$dira=pi-data$diry*pi/180 + dira_offset
			data$dire=pi/2-data$dirx*pi/180
			}
		else if(sonR_implemented(data, c("SBE","MBE"))){
			# Corresponding to the method used in EKRaw2TSD.m:
			data$dira=data$dirx*pi/180 + dira_offset
			data$dire=pi-data$diry*pi/180
			# Turn the beams that have elevation angle greater than pi by the angle pi in azimuth angle, and then mirror the elevation angles around pi:
			toboturned=data$dire>pi
			data$dira[toboturned]=(pi+data$dira[toboturned]) %% (2*pi)
			data$dire[toboturned]=pi-(data$dire[toboturned]-pi)
			}
		else if(sonR_implemented(data, "OFS")){
			# Corresponding to the method used in EKRaw2TSD.m:
			# data$dira=((data$diry)*pi/180) %% (2*pi)
			data$dira=pi - ((data$diry)*pi/180) %% (2*pi) + dira_offset
			data$dire=pi/2 + data$dirx*pi/180
			}
		else{
			stop(paste("Unknown acoustical instrument specified by 'esnm' (",data$esnm[1],")",sep=""))
			}
		}
	data
	##################################################
	##################################################
	}
