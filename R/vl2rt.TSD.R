#*********************************************
#*********************************************
#' Transforms from velocity data to rotation data, given as counter clockwise around z-rotation and counter clockwise around x-rotation.
#'
#' @param data  is a list containing the variables 'vlx*', 'vly*' and 'vlz*'.
#' @param var  is a vector of the elements to return.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD strff
#'
#' @export
#' @rdname vl2rt.TSD
#'
vl2rt.TSD<-function(data,var=c("rtz","rtx"),fresh=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-12-05 - Clean version.
	########### DESCRIPTION: ###########
	# Transforms from velocity data to rotation data, given as counter clockwise around z-rotation and counter clockwise around x-rotation.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing the variables 'vlx*', 'vly*' and 'vlz*'.
	# ---var--- is a vector of the elements to return.
	

	##################################################
	##################################################
	##### Preparation and execution #####
	if((fresh || length(data$rtzf)!=length(data$vlxf)) && strff("rtz",var)){
		if(is.list(data$vlyf)){
			data$rtzf=data$vlyf
			for(i in seq_along(data$vlyf)){
				data$rtzf[[i]] = (atan2(data$vlyf[[i]],data$vlxf[[i]])-pi/2) %% (2*pi)
				}
			}
		else if(length(data$vlyf)>0){
			data$rtzf = (atan2(data$vlyf,data$vlxf)-pi/2) %% (2*pi)
			}
		else{
			stop("Rotation angle around the z axis unavailable in 'x'")
			}
		}
	if((fresh || length(data$rtxf)!=length(data$vlxf)) && strff("rtx",var)){
		if(is.list(data$vlyf)){
			data$rtxf=data$vlyf
			for(i in seq_along(data$vlyf)){
				disthorizontal = sqrt(data$vlxf[[i]]^2+data$vlyf[[i]]^2+data$vlzf[[i]]^2)
				data$rtxf[[i]] = pi/2 - acos(data$vlzf[[i]]/disthorizontal)
				# Default missing angles to 0
				data$rtxf[[i]][is.na(data$rtxf[[i]])] = 0
				}
			}
		else if(length(data$vlyf)>0){
			disthorizontal = sqrt(data$vlxf^2+data$vlyf^2+data$vlzf^2)
			data$rtxf = pi/2 - acos(data$vlzf/disthorizontal)
			# Default missing angles to 0
			data$rtxf[is.na(data$rtxf)] = 0
			}
		else{
			stop("Rotation angle around the x axis unavailable in 'x'")
			}
		}


	##### Output #####
	data
	##################################################
	##################################################
	}
