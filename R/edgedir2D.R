#*********************************************
#*********************************************
#' Calculates the edge points of an acoustic device based on the CTD-data, using edge points along the beam maxima as basis, and extracting the vertical dimension (elevation angle) from the vertical angular distance between beam maxima. !!!Only one time step is treated!!!!
#'
#' @param data  is the list of TSD inputs as returned from read.TSD (must contain ""dira" and "dire", and CTD-data).
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
#' @importFrom TSD ang2rot car2sph NAs strff zeros
#'
#' @export
#' @rdname edgedir2D
#'
edgedir2D<-function(data, seabed=-12000, rot=1, compensation=c("pitch","roll"), ideal=TRUE, stretch=1, fanWidth="b2", ...){
	#edgedir2D<-function(data, esnm=NULL, seabed=-12000, rot=1, compensation=c("pitch","roll"), ideal=TRUE, stretch=1, fanWidth="b2", ...){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2015-03-26 - Clean version.
	# Last: 2016-10-28 - Fixed bug with ME70. Before only the beam maxima were rotated to a fan facing in the direction of the vessel, but also the beam direction angles need to be treated in the same way.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Extract the name of the acoustic instrument:
	#if(is.null(esnm)){
	#	esnm=data$esnm
	#	}
	#else if(!identical(tolower(data$esnm),tolower(esnm))){
	#	warning("The value of 'esnm' given in 'data' ignored")
	#	}
	# If 'rect' is not in 'data', issue a warning:
	if(length(data$rect)==0){
		warning("The list 'data' must contain the logical vector 'rect' defining which beams are rectangular (FALSE for circular/elliptical beams)")
		}
	
	##### Execution #####
	# All one-fan multibeam sonars and echosounders are treated in such a way that the widths of the voxels are determined from the distance between beams along the fan and the height of the voxels are set to the mean of the heights of each beam corresponding to the equivalent beam angle. The ME70 multibeam echosounder points downwards and but must undergo the following procedure to ensure that beam refraction and the correct volume is included: (1) Rotate the beam directions to aim along the sea surface instead of vertical. (2) Expand the beam directions to the edges instead of mid points of voxels (3) Rotate back to facing vertical downwards. (4) Run soundbeam.TSD on the edges of the beams. (5) Rotate to the sea surface again. For all other one-fan instruments steps (2) and (4) are done:
	
	### THIS METHOD WAS UNSUCCESSFULLY REPLACED BY A METHOD THAT TRACED THE SOUND BEAMS ALONG THE EDGES, WHICH WAS INTENDED TO EASE THE CODE AND UNDERSTANDING OF THE CODE. THE REASON WHY THIS WAS UNSUCCESSFUL WAS THAT FOR THE ME70 SONAR THE RANGE TO THE VOXELS IS USED TO CREATE WIDER VOXELS ALONG THE DIRECTION OF MOTION OF THE VESSEL WHICH STRETCH != 1, WHICH SMOOTHS THE PPLOT. THUS, THE RANGE TO THE VOXELS NEEDS TO BE CALCULATED FIRST, AS DONE HERE, AND THEN THE VOXELS ARE SPANNED FROM THE BEAM MAXIMA AND OUT TO IN THE AZIMUTH AND THE ELEVATION DIRECTION. USING THE OLD METHOD HERE:
	
	### (1) ###
	out_r = car2sph( soundbeam.TSD(data, ind=data$rect, cs="v", seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, rpos="edge", drop.out=FALSE, plot=FALSE), list.out=FALSE, nonneg=TRUE)
	
	### (2) ###
	# Rotate the beams to point in the direction of the vessel if vertical is included in the fan. This includes first rotating the coordinate system of the vessel around the z axis so that the x axis coincides with the azimuth angle, and then rotating around the x axis by negative 90 degrees, so point the y axis in the direciton of the fan (downwards):
	#out_r = rotate3D(out_r, by=data$offset$by, ang=data$offset$ang, sph.in=TRUE, sph.out=TRUE, list.out=FALSE)
	vertical = any(data$dire > 170*pi/180)
	if(vertical){
		out_r = rotate3D(out_r, by=data$offset$by, ang=data$offset$ang, sph.in=TRUE, sph.out=TRUE, list.out=FALSE)
		}
	
	# Set dimensions only to discard the closest samples in the following code:
	D1 = max(data$lenb[data$rect])+1
	#D3 = length(unique(data$dire[data$rect]))
	D2 = sum(data$rect)
	dim(out_r) = c(D1, D2, 3)
	# Convert from angles to rotation angles moving beyond the range 0, 2*pi:
	out_r[-1,,2] = ang2rot(out_r[-1,,2,drop=FALSE])
	out_r[-1,,3] = ang2rot(out_r[-1,,3,drop=FALSE])
	
	
	### (3) ###
	# Get the lower and upper edge r-values of the voxels:
	r1 = c(out_r[-1,,1])
	r0 = c(out_r[-D1,,1])
	
	### (4) ###
	##dataRot = rotate3D(cbind(1, data$dira, data$dire), by=data$offset$by, ang=data$offset$ang, sph.in=TRUE, sph.out=TRUE, list.out=FALSE)
	##data$dire = dataRot[,2]
	##data$dira = dataRot[,3]
	
	
	
	### For the treatment of angles, consider the far edge of each voxel, by out_r[-1, ...]. This also excludes possible zero values at the first range of the sonar volume: ###
	# Get the lower and upper edge theta-values of the voxels, using the difference in azimuth angle between voxels:
	halfdiff = zeros(D1-1, D2)
	halfdiff[,-1] = (out_r[-1,-1,2,drop=FALSE] - out_r[-1,-D2,2,drop=FALSE]) / 2
	halfdiff[,1] = halfdiff[,2]
	# Expand to one half theta on each side:
	theta1 = out_r[-1,,2] - halfdiff
	theta0 = out_r[-1,,2] + halfdiff
	# Convert back to 0, 2*pi:
	#theta1 = theta1 %% (2*pi)
	#theta0 = theta0 %% (2*pi)
	
	
	### (5) ###
	# Get the height of the fan:
	data <- getFanWidth(data=data, fanWidth=fanWidth, ind=data$rect, stretch=stretch, range=r1)
	
	# Expand to one half phi on each side:
	phi1 = out_r[-1,,3] - rep(data$diffdire/2, each=D1-1)
	phi0 = out_r[-1,,3] + rep(data$diffdire/2, each=D1-1)
	#phi1 = c(out_r[-1,,3]) - c(data$diffdire/2)
	#phi0 = c(out_r[-1,,3]) + c(data$diffdire/2)
	# Convert back to 0, 2*pi:
	phi1 = phi1 %% pi
	phi0 = phi0 %% pi
	
	# Combine the above vectors:
	out_r = cbind(r0=c(r0), r1=c(r1), theta0=c(theta0), theta1=c(theta1), phi0=c(phi0), phi1=c(phi1))
	
	
	##### Output #####
	# Return a list of the edgepoints of the rectangular voxels:
	return(list(rect=out_r))
	##################################################
	##################################################
	}
