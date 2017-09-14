#*********************************************
#*********************************************
#' Calculates the volumes of the voxels of an omnidirectional fan of a multibeam sonar. Only the rectangular beams are treated!!! See read.event() for importing volumes of voxels.
#'
#' @param data  is the list of the inputs used in the function as returned from read.TSD() (must include "").
#' @param fanWidth  has a number of possible values (abbreviations allowed, so that the values "e" and "b" can be used):
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#'
#' @export
#' @rdname volx2D_oneping
#'
volx2D_oneping<-function(data, fanWidth="b2"){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-07-24 - Clean version.
	# Last: 2015-11-24 - Return list.
	########### DESCRIPTION: ###########
	# Calculates the volumes of the voxels of an omnidirectional fan of a multibeam sonar. See read.event() for importing volumes of voxels.
	########## DEPENDENCIES: ###########
	# zeros() 
	############ VARIABLES: ############
	# ---data--- is the list of the inputs used in the function as returned from read.TSD() (must include "").
	# ---fanWidth--- has a number of possible values (abbreviations allowed, so that the values "e" and "b" can be used):
	#	1) "equivalentBeamangle", causing the widths of the voxels across the fan to be calculated independently, using the widhts along the fan (calculated by the angle between the beam maxima) and the equivalent beam angle.
	#	2) "beamwidth" using the across beam beam width as the across beam voxel width.
	#	3) numeric, setting the across fan voxel width directly.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# The volume of the spherical segments, calculated by the following argument:
	# Volume of an open spherical sector (see http://mathworld.wolfram.com/SphericalSector.html), is 2 pi r^3(cos(phi1) − cos(phi2))/3, where 'r' is the radial distance to the edge of the sphere, and 'theta' and 'phi' are the azimuth and elevation angles, respectively. Subtracting for two radial distances and subsetting the resulting strip by theta relative to 2*pi gives: 
	# (theta2−theta1) (cos(phi1) −cos(phi2))(r2^3 −r1^3) / 3:
	
	# If average speed of sound is not given, it is set to 1500 with a warning.
	if(is.null(data$asps) && is.null(data$asps)){
		warning("Average speed of sound not given, and is set to 1500 m/s")
		data$asps=1500
		}
	
	# Get the equivalent beam angle in radians. The equivalent beam angle is defined in SIMRAD systems as 10*log10(psi), where 'psi' is the convetional equivalent beam angle (see "Simrad, EK500 Power to Sv and TS.pdf"):
	#data = getVoxelWidth2D(data, fanWidth=fanWidth)
	data = getFanWidth(data, fanWidth=fanWidth, ind=data$rect, stretch=1)
	# Calculate the radii of the ranges used in the calculation 
	data$lenb = max(data$lenb)
	#dr=data$asps*data$sint/2
	#ranges = c(0, 1:data$lenb*dr-dr/2)
	data$ranges = soundbeam_range(data, pos="edge")
	# Partition into open spherical sectors:
	data$volx = outer( diff(data$ranges^3), abs(cos(data$direExpanded[,1]) - cos(data$direExpanded[,2])) * data$diffdira/3 )
	

	##### Output #####
	data[c("volx", "ranges", "diffdira", "diffdire", "direExpanded")]
	##################################################
	##################################################
	}
