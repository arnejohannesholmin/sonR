#*********************************************
#*********************************************
#' Calculates the volumes of the voxels of the EK60 multifrequency echosounder given the CTD-data and beams-dataof the particular data set 'dataset'. See read.event() for importing volumes of voxels.
#'
#' @param data  is the list of the inputs used in the function as returned from read.TSD() (must include "").
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom TSD ind.expand
#'
#' @export
#' @rdname volx1D_oneping
#'
volx1D_oneping<-function(data){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-07-24 - Clean version.
	# Update: 2010-02-09 - Added support for TSD input and adapted to the array structure having the radial part along the first dimension and the beams along the second (lenb,numb).
	# Last: 2015-11-24 - Return list.
	########### DESCRIPTION: ###########
	# Calculates the volumes of the voxels of the EK60 multifrequency echosounder given the CTD-data and beams-dataof the particular data set 'dataset'. See read.event() for importing volumes of voxels.
	########## DEPENDENCIES: ###########
	# zeros()
	############ VARIABLES: ############
	# ---data--- is the list of the inputs used in the function as returned from read.TSD() (must include "").
	
	
	##################################################
	##################################################
	##### Preparation #####
	# The volume of circular voxels is 2/3*pi*(1-cos(eqba/2)*(r_(i) - r_(i-1))), found from "http://mathworld.wolfram.com/SphericalCone.html":
	# If average speed of sound is not given, it is set to 1500 with a warning.
	if(is.null(data$asps) && is.null(data$asps)){
		warning("Average speed of sound not given, and is set to 1500 m/s")
		data$asps=1500
		}
	
	# Get the equivalent beam angle in radians. The equivalent beam angle is defined in SIMRAD systems as 10*log10(psi), where 'psi' is the convetional equivalent beam angle (see "Simrad, EK500 Power to Sv and TS.pdf"):
	eqba=10^(data$eqba/10)
	data$lenb=max(data$lenb)
	data$numb=max(data$numb)
	#dr=data$asps*data$sint/2
	
	# Apply the subset of the beams given by 'ind':
	ind = ind.expand(!data$rect, data$numb, drop=TRUE)
	data$eqba=10^(data$eqba[ind]/10)
	
	
	##### Execution #####
	#ranges=c(0,1:data$lenb*dr-dr/2)
	data$ranges = soundbeam_range(data, pos="edge")
	# Transform to half opening angle of the spherical sector, using the expression found on "http://en.wikipedia.org/wiki/Steradian"
	data$theta=acos(1-data$eqba/2/pi)
	# The volume of a round voxel is 2/3*pi*(1-cos(eqba/2)*(r_(i) - r_(i-1))), found from "http://mathworld.wolfram.com/SphericalCone.html":
	data$volx=2/3*pi * outer( diff(data$ranges^3) , (1-cos(data$theta)) )
	
	
	##### Output #####
	data[c("volx", "ranges", "theta")]
	##################################################
	##################################################
	}
