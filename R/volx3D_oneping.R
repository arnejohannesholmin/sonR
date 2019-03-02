#*********************************************
#*********************************************
#' Calculates the volumes of the voxels of the MS70 multibeam sonar given the CTD-data and beams-data. See read.event() for importing volumes of voxels.
#'
#' @param data  is the list of the inputs used in the function as returned from read.TSD() (must include "").
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#'
#' @export
#' @rdname volx3D_oneping
#'
volx3D_oneping<-function(data){
		
	############### LOG: ###############
	# Start: 2009-07-24 - Clean version.
	# Update: 2010-02-09 - Added support for TSD input and adapted to the array structure having the radial part along the first dimension and the beams along the second (lenb,numb).
	# Last: 2015-11-24 - Return list.
	
	##### Preparation #####
	# The volume of the spherical segments, calculated by the following argument:
	# Volume of an open spherical sector (see http://mathworld.wolfram.com/SphericalSector.html), is 2 pi r^3(cos(phi1) - cos(phi2))/3, where 'r' is the radial distance to the edge of the sphere, and 'theta' and 'phi' are the azimuth and elevation angles, respectively. Subtracting for two radial distances and subsetting the resulting strip by theta relative to 2*pi gives: 
	# (theta2 - theta1) (cos(phi1) -cos(phi2))(r2^3 - r1^3) / 3:
	
	# If average speed of sound is not given, it is set to 1500 with a warning.
	if(is.null(data$asps) && is.null(data$asps)){
		warning("Average speed of sound not given, and is set to 1500 m/s")
		data$asps=1500
		}
	
	lenb=max(data$lenb)
	numb=max(data$numb)
	#dr=data$asps[1]*data$sint[1]/2
	
	# The unique azimuth angles:
	data$udira=sort(unique(data$dira))
	# The unique elevation angles and indexes for which horizontal fan each beam belongs to:
	data$udire=rev(sort(unique(data$dire)))
	inde=match(data$dire, data$udire)
	# Angular partitioning:
	da=data$udira[2]-data$udira[1]
	de=data$udire[2]-data$udire[1]
	# Vector of elevation angles for the partitioning surfaces between horizontal fans:
	data$udire=c(data$udire[1]-de/2, data$udire+de/2)
	
	
	##### Execution #####
	# The output:
	data$ranges=soundbeam_range(data, pos="edge")
	
	#ranges=c(0,1:lenb*dr-dr/2)
	data$volx=outer( diff(data$ranges^3), rep(diff(cos(data$udire)) * da/3, each=length(data$udira)) )
	
	
	##### Output #####
	# Horizontal partitioning:
	data[c("volx", "ranges", "udira", "udire")]
}
