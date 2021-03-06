#*********************************************
#*********************************************
#' Calculates the horizontal area of the voxels of the SH80 multibeam sonar given the CTD-data and beams-data of the particular data set 'dataset'. See read.event() for importing volumes of voxels.
#'
#' @param data  is the list of the inputs used in the function as returned from read.TSD() (must include "").
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom utils tail
#'
#' @export
#' @rdname harx.SH80_oneping
#'
harx.SH80_oneping<-function(data){
		
	# The horizontal area of the spherical segments, calculated by the following argument:
	# The difference in area of a circle at two distances cos(r1) and cos(r2) (where 'r1' and 'r2' are taken at the intersection between the acoustic axis of the beams and the spheres at radii 'r1' and 'r2') multiplied by the fractional angle of the beams compared to the whole circle:
	# If average speed of sound is not given, it is set to 1500 with a warning.
	if(is.null(data$asps) && is.null(data$asps)){
		warning("Average speed of sound not given, and is set to 1500 m/s")
		data$asps = 1500
		}
	lenb = max(data$lenb)
	
	
	##### Execution #####
	# Get the azimuth angle differences:
	# If only one beams is present in the data, assume that it spans a full circle in azimuth angle (in the x-y-plane):
	if(length(data$dira)==1){
		diffdira <- 2 * pi
	}
	else{
		diffdira = diff(data$dira)
		diffdira = c(diffdira,tail(diffdira,1))
	}
	
	# Calculate the radii of the spheres used in the calculation:
	spheres = soundbeam_range(data, pos="edge") * cos(data$dire[1]-pi/2)
	#dr=data$asps[1]*data$sint[1]/2
	#spheres=c(0,1:lenb*dr-dr/2) * cos(data$dire[1]-pi/2)
	# Partition into radial fractions of circle segments:
	out = outer(diff(spheres^2), abs(diffdira)/2)
	
	
	return(out)
}
