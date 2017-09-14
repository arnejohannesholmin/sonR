#*********************************************
#*********************************************
#' Transforms from 'sgsc' to positions, in a similar way as used in pplot3d.TSD(). Used in read.event() for clustering 'sgsc'.
#'
#' @param data  is the list of TSD inputs as returned from read.TSD(var=c("vbsc","voxels","vessel")). Should only contain one ping/timestep.
#' @param N  is the approximate number of points plotted.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param ideal  is TRUE to represent the simple case where the speed of sound 'data$asps' is invariant of depth.
#' @param allert  is the limit of the number of points inducing a warning asking the user to continue. Default is NULL, indicating no check for the number of points to be generated.
#' @param seabed  is the depth of the seabed, at which the beams are reflected when calculating the midpoints of voxels.
#' @param rot  see soundbeam.TSD().
#' @param compensation  is a vector of string giving which rotation values that are compensated for in the sonar. Only c("pitch","roll") is available for the current version. Used in soundbeam.TSD.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#'
#'
#'
#'
#' @export
#' @rdname read.event_sgsc2pos.TSD
#'
read.event_sgsc2pos.TSD<-function(data,N=1e2,esnm="MS70",ideal=TRUE,allert=NULL,seabed=-12000,rot=1,compensation=c("pitch","roll"),...){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-08-28 - Clean version.
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(!is.null(data$esnm)){
		esnm=data$esnm
		}
	# Assure that the required variables are present:
	requiredVariables1=c("dira","dire","lenb","esnm","numb","eqba")
	requiredVariables2=c("sgsc","sgs0")
	if(!(all(requiredVariables1 %in% names(data)) & any(requiredVariables2 %in% names(data)))){
		stop("Segmentation data 'sgsc'/'sgs0', must be present in 'data'")
		}
	
	# Get the edges of the voxels (rectangular or circular) before subsetting the data:
	data$psxv=0
	data$psyv=0
	data$pszv=0
	data$rtzv=0
	data$psze=0
	rthetaphi=voxel.edge.TSD(data, t=1, esnm=esnm, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal)[[1]]
	# Select only the voxels in 'sgsc':
	rthetaphi=rthetaphi[data$sgsc,]
	
	
	##### Execution #####
	# Select the voxels in which more than one points is to be generated:
	n=rep(N,length(data$sgsc))
	if(length(allert)>0){
		N=sum(n,na.rm=TRUE)
		if(N>allert){
			ask=identical("y", tolower(readline(paste(N," positions will be generated. Continue? ('y')"))))
			if(!ask){
				return()
				}
			}
		}
		
	# Regenerate the school by drawing spherically uniformly distributed points within each voxel in the three cases: (1) all points, (2) weak voxels, (3) voxels with more than nlim points:
	regen=runif.sph(n, r=rthetaphi[,1:2], theta=rthetaphi[,3:4], phi=rthetaphi[,5:6], sph.out=FALSE, offset=get.specs.esnm(data)$offset)
	
	
	##### Output #####
	regen
	##################################################
	##################################################
	}
