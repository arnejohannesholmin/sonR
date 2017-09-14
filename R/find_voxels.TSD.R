#*********************************************
#*********************************************
#' Finds the voxels in which the positions in 'pos' are located. NAs returned for voxels outside of the sonar volume.
#'
#' @param pos  is a three column matrix of the r, theta, and phi spherical position of the points to locate to voxels.
#' @param data  is the list of beam configuration data.
#' @param arr.ind  is TRUE to return array indices, in which voxels outside of the sonar volume are indicated by zeros in one oner more of the columns.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom TSD arr.ind2ind
#'
#' @export
#' @rdname find_voxels.TSD
#'
find_voxels.TSD<-function(pos, data, arr.ind=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-03-23 - Clean version.
	# Last: 2014-08-06 - Added 'arr.ind'.
	########### DESCRIPTION: ###########
	# Finds the voxels in which the positions in 'pos' are located. NAs returned for voxels outside of the sonar volume.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---pos--- is a three column matrix of the r, theta, and phi spherical position of the points to locate to voxels.
	# ---data--- is the list of beam configuration data.
	# ---arr.ind--- is TRUE to return array indices, in which voxels outside of the sonar volume are indicated by zeros in one oner more of the columns.
	
	
	##################################################
	##################################################
	##### Preparation #####
	dims <- c(max(data$lenb), length(unique(data$dira)), length(unique(data$dire)))
	
	
	##### Execution and output #####
	#data$dr <- data$sint*data$asps/2
	# Indexes locating the fish at the correct radial layer, assuming constant speed of sound:
	#r <- findInterval(pos[,1]/data$dr,0:data$lenb[1])
	r <- findInterval(pos[,1], soundbeam_range(data, pos="edge"))
	
	data$dira <- unique(data$dira)
	data$dtheta <- diff(data$dira[1:2])
	data$dira <- c(data$dira[1] - data$dtheta, data$dira) + data$dtheta/2
	reva <- FALSE
	if(data$dira[2]<data$dira[1]){
		reva <- TRUE
		data$dira <- rev(data$dira)
		}
	data$dire <- unique(data$dire)
	data$dphi <- diff(data$dire[1:2])
	data$dire <- c(data$dire[1] - data$dphi, data$dire) + data$dphi/2
	reve <- FALSE
	if(data$dire[2]<data$dire[1]){
		reve <- TRUE
		data$dire <- rev(data$dire)
		}
	
	theta <- findInterval(pos[,2]%%(2*pi), data$dira)
	if(reva){
		theta <- length(data$dira) - theta
		}
	phi <- findInterval(pos[,3], data$dire)
	if(reve){
		phi <- length(data$dire) - phi
		}
	
	# Convert to voxel indexes:
	if(!arr.ind){
		arr.ind2ind(cbind(r, theta, phi), dims)
		}
	else{
		cbind(r, theta, phi)
		}
	##################################################
	##################################################
	}
