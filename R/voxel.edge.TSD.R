#*********************************************
#*********************************************
#' Calculates a matrix of edge points of voxels of an underwater acoustic system. The matrix contains r0 and r1, the closer and father range to each voxel; theta0 and theta1, the azimuth angle of the left and right dge of each voxel; and phi0 and phi1, the elevation angle of the lower and upper edge of each voxel, as seen from the vessel for a multi-beam sonar such as the Simrad MS70. !!!Only one time step is treated!!!!
#'
#' @param data  is the list of inputs variables as returned from read.TSD (including ).
#' @param t  is a single integer giving the time step for which the edge matrix in the coordinate system of the vessel is returned.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param seabed  is the z-coordinate of the sea bed, usually provided by echo sounder data. Soundbeams reflected from the sea bed or the surface are reflected at the incidence angle.
#' @param rot  is 1 if simple rotation using cosine and sine is to be used, and 2 if the function rotate() is to be used in the rotation. Times for the different methods (tested on MacBook Pro dual 2.8 GHz, 2010-02-09, with other applications running):
#' @param compensation  is a vector of string giving which rotation values that are compensated for in the sonar. Only c("pitch","roll") is available for the current version. Used in soundbeam.TSD.
#' @param ideal  is TRUE to represent the simple case where the speed of sound 'data$asps' is invariant of depth.
#' @param stretch  is used to stretch the voxels of the ME70 multibeam echosounder in the direction of motion, so that space in between voxels is smoothed out.
#' @param fanWidth  has a number of possible values: (1) "b1": one way beam width. (2) "b2": two way beam width. (3) "fe": beams modeled by rectangular cones with width withing the fan given by the inter-beam angle, and calculated using the equivalent beam angle. This normally causes larged fan width due to overlap between beams.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname voxel.edge.TSD
#'
voxel.edge.TSD <- function(data, t=1, esnm=NULL, seabed=-12000, rot=1, compensation=c("pitch", "roll"), ideal=TRUE, stretch=1, fanWidth="b2", ...){
	
	############### LOG: ###############
	# Start: 2011-09-25 - Clean version.
		
	##### Preparation #####
	# Function for selecting one ping of multiple pings beam configuration data (particularly for fishery sonars):
	getPingBeams = function(x, p){
		c(lapply(x[c("asps","sint","numb")], function(xx) xx[p]), lapply(x[c("lenb","dira","dire","eqba")], function(xx) xx[,p]))
		}
	# Extract the name of the acoustic instrument:
	if(!is.null(data$esnm[1])){
		esnm=data$esnm[1]
		}
	data = get.specs.esnm(data)
	
	
	##### Execution and output #####
	# Get the edges of the voxels, rectangular and/or circular:
	if(sonR_implemented(esnm, "SBE")){
		edgedir1D(data, esnm=esnm, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, ...)
		}
	else if(sonR_implemented(esnm, "MBE")){
		c(	edgedir2D(data, esnm=esnm, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, stretch=stretch, fanWidth=fanWidth, ...), 
			edgedir1D(data, esnm=esnm, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, stretch=stretch, fanWidth=fanWidth, ...))
		# Discard the split beams:
		###edgedir2D(data, esnm=esnm, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, stretch=stretch, fanWidth=fanWidth, ...)
		}
	else if(sonR_implemented(esnm, "MBS")){
		edgedir3D(data, esnm=esnm, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, ...)
		}
	else if(sonR_implemented(esnm, "OFS")){
		# Here only one time step is treated, since the function that used voxel.edge.TSD(), pplot3d_sv2pos.TSD(), only considers one time step. Thus the beam configuration is extracted for the current time step in that function:
		edgedir2D(data, esnm=esnm, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, stretch=stretch, fanWidth=fanWidth, ...)
		}
	else{
		warning(paste("Unsupported echo sounder ",esnm," for voxel data",sep=""))
		}
	}
