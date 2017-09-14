#*********************************************
#*********************************************
#' Calculates edge values for voxels in a 2-D sonar fan. If a vertically oriented fan is 
#'
#' @param data  is the list of the inputs used in the function as returned from read.TSD() (must include "").
#' @param fanWidth  has a number of possible values: (1) "b1": one way beam width. (2) "b2": two way beam width. (3) "fe": beams modeled by rectangular cones with width withing the fan given by the inter-beam angle, and calculated using the equivalent beam angle. This normally causes larged fan width due to overlap between beams.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ind.expand NAs strff zeros
#'
#' @export
#' @rdname getVoxelWidth2D
#'
getVoxelWidth2D<-function(data, fanWidth="b2"){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-07-24 - Clean version.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Use only rectangular beams:
	#numb = max(data$numb)
	#ind = ind.expand(data$rect, numb, drop=TRUE)
	#beamsVar <- sapply(data, length) == length(data$dira)
	#data[beamsVar] = lapply(data[beamsVar], "[", ind)
	#data$dira = data$dira[ind]
	#data$dire = data$dire[ind]
	#data$eqba = data$eqba[ind]
	#eqba = 10^(data$eqba[ind]/10)
	
	
	##### Execution #####
	# Identify whether the beams cover the vertical direction, by assuming that i one beam has elevation angle larger than 180-10 degrees, then the fan covers "vertical":
	##vertical = any(data$dire > 170*pi/180)
	# Rotate the beams to point in the direction of the vessel if vertical is included in the fan. This includes first rotating the coordinate system of the vessel around the z axis so that the x axis coincides with the azimuth angle, and then rotating around the x axis by negative 90 degrees, so point the y axis in the direciton of the fan (downwards):
	##if(vertical){
		#dirRotated = rotate3D(cbind(1, data$dira, data$dire), by=paste("z",data$offset$by), ang=c(data$dira[1],data$offset$ang), sph.in=TRUE, sph.out=TRUE, list.out=FALSE)
		##dirRotated = rotate3D(cbind(1, data$dira, data$dire), by=data$offset$by, ang=data$offset$ang, sph.in=TRUE, sph.out=TRUE, list.out=FALSE)
		##data$dira = dirRotated[,2]
		##data$dire = dirRotated[,3]
		##}
	
	
	# Get the angles of the rectangular voxels based on the equivalent beam angles and the directions of the beam maxima:
	#data$diffdira = abs(diff(data$dira))
	#data$diffdira = c(data$diffdira, tail(data$diffdira,1))
	
	data <- getFanWidth(data=data, fanWidth=fanWidth, ind=data$rect, stretch=1)
		
	# Add half the voxel width on each side of the elevation direciton angles:
	#data$direExpanded = cbind(data$dire-data$diffdire/2 , data$dire+data$diffdire/2)
	
	
	##### Output #####
	# Return the voxel widths:
	data
	##################################################
	##################################################
	}
