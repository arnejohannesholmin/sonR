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
#' @param ... Used for robutness.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD car2sph zeros
#'
#' @export
#' @rdname edgedir3D
#'
edgedir3D<-function(data, esnm=NULL, seabed=-12000, rot=1, compensation=c("pitch","roll"), ideal=TRUE, ...){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2015-03-26 - Clean version.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Extract the name of the acoustic instrument:
	if(is.null(esnm)){
		esnm=data$esnm
		}
	else if(!identical(tolower(data$esnm),tolower(esnm))){
		warning("The value of 'esnm' given in 'data' ignored")
		}
	# If 'rect' is not in 'data', issue a warning:
	if(length(data$rect)==0){
		warning("The list 'data' must contain the logical vector 'rect' defining which beams are rectangular (FALSE for circular/elliptical beams)")
		}
	
	
	##### Execution #####
	# Three dimensional sonar requires a spherically data$rectangular grid of equal voxels sizes in the spherical coordinate system, implying decreasing voxel sizes in the cartesian coordinate system as the elevation angle approaches the poles:
	out_r = car2sph( soundbeam.TSD(data, ind=data$rect, cs="v", seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, rpos="edge", drop.out=FALSE, plot=FALSE), list.out=FALSE, nonneg=TRUE)
	D1 = max(data$lenb[data$rect])+1
	D3 = length(unique(data$dire[data$rect]))
	D2 = length(data$dire[data$rect]) / D3
	dim(out_r) = c(D1, D2, D3, 3)
	
	# Get the lower and upper edge r-values of the voxels:
	r1 = c(out_r[-1,,,1])
	r0 = c(out_r[-D1,,,1])
	
	
	### For the treatment of angles, consider the far edge of each voxel, by out_r[-1, ...]: ###
	# Get the lower and upper edge theta-values of the voxels, using the difference in azimuth angle between voxels:
	halfdiff = zeros(D1-1, D2, D3)
	halfdiff[,-1,] = (out_r[-1,-1,,2] - out_r[-1,-D2,,2]) / 2
	halfdiff[,1,] = halfdiff[,2,]
	# Expand to one half theta on each side:
	theta1 = out_r[-1,,,2] - halfdiff
	theta0 = out_r[-1,,,2] + halfdiff
	
	
	# Get the lower and upper edge theta-values of the voxels, using the difference in azimuth angle between voxels:
	halfdiff = zeros(D1-1, D2, D3)
	halfdiff[,,-1] = (out_r[-1,,-1,3] - out_r[-1,,-D3,3]) / 2
	halfdiff[,,1] = halfdiff[,,2]
	# Expand to one half theta on each side:
	phi1 = out_r[-1,,,3] - halfdiff
	phi0 = out_r[-1,,,3] + halfdiff
	
	
	# Combine the above vectors:
	out_r=cbind(r0=c(r0), r1=c(r1), theta0=c(theta0), theta1=c(theta1), phi0=c(phi0), phi1=c(phi1))
	
	
	##### Output #####
	# Return a list of the edgepoints of the rectangular voxels:
	return(list(rect=out_r))
	##################################################
	##################################################
	}
