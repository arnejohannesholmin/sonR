#*********************************************
#*********************************************
#' Applies a Gaussian 1-D, 2-D, or 3-D kernel on the 'data' given as a list with names according to the TSD-convension, and the bandwidths given in the vector 'h':
#'
#' @param data  is the list of data with names according to the TSD-convension.
#' @param h  is the three element vector of bandwidths.
#' @param w  is the boundary of the kernel, outside which the it is 0.
#' @param sim  is a TRUE if smoothing should be done only along the first dimensions, simultaneously over the stages of the last dimension. If 'sim' is an integer larger than 1, the positions 'coords' are used 'sim' times, and the data 'x' should have length 'sim' times the length of one coordinate of 'coords'.
#' @param ind  is a list of indexes, as given to subset_TSD(), used to select the subset over which the estimation of high intensity noise is done. Defaulted to exclude the first 100 voxels along each beam.
#' @param drop  is TRUE if the data should be dropped of unused dimensions prior to smoothing. This can be useful for smoothing data which .
#' @param var  is a string specifying the variable to smooth.
#' @param na.rm  is single integer representing the dimension along which NAs are discarded from the smoothing in the case that sim==1. For example, if na.rm=2 and the dimension of 'x' is [5,12,7], and x[3,2:4,5]=NA, then all data x[,2:4,] will be excluded from the smoothing and set to NA. If na.rm=FALSE, no NAs should be contained in the data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom utils tail
#'
#' @export
#' @rdname ksmooth.SG.TSD
#'
ksmooth.SG.TSD <- function(data, h=1, w=h*3, sim=FALSE, ind=list(), drop=FALSE, var="vbsc", na.rm=1){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-08-23 - Collected several functions into one (like kernSmooth3dGaussMultipleNoNA.TSD() and kernSmoothGauss()).
	########### DESCRIPTION: ###########
	# Applies a Gaussian 1-D, 2-D, or 3-D kernel on the 'data' given as a list with names according to the TSD-convension, and the bandwidths given in the vector 'h':
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is the list of data with names according to the TSD-convension.
	# ---h--- is the three element vector of bandwidths.
	# ---w--- is the boundary of the kernel, outside which the it is 0.
	# ---sim--- is a TRUE if smoothing should be done only along the first dimensions, simultaneously over the stages of the last dimension. If 'sim' is an integer larger than 1, the positions 'coords' are used 'sim' times, and the data 'x' should have length 'sim' times the length of one coordinate of 'coords'.
	# ---ind--- is a list of indexes, as given to subset_TSD(), used to select the subset over which the estimation of high intensity noise is done. Defaulted to exclude the first 100 voxels along each beam.
	# ---drop--- is TRUE if the data should be dropped of unused dimensions prior to smoothing. This can be useful for smoothing data which .
	# ---var--- is a string specifying the variable to smooth.
	# ---na.rm--- is single integer representing the dimension along which NAs are discarded from the smoothing in the case that sim==1. For example, if na.rm=2 and the dimension of 'x' is [5,12,7], and x[3,2:4,5]=NA, then all data x[,2:4,] will be excluded from the smoothing and set to NA. If na.rm=FALSE, no NAs should be contained in the data.
		

	##################################################
	##################################################
	##### Preparation #####
	# Issue a warning if the specified variable to smooth has zero length:
	if(length(data[[var]])==0){
		warning(paste("The specified variable \"", var, "\" is not present in the data", sep=""))
		}
	# Check that beam configuration data are present:
	if(any(length(data$lenb)==0, length(data$numb)==0, length(data$freq)==0)){
		stop("Beam configuration data 'lenb', 'numb' and 'freq' must be present in the data")
		}
	
	# If the length of the variable specified by 'var' is a multiple > 1 times the length of the positions data, 'sim' is set to TRUE, and the same positions are used for all time steps:
	lx = length(data[[var]])
	lcoords = length(data$psxx)
	if(sim==0 && ((lx/lcoords) %% 1)==0){
		if(lx>lcoords){
			warning("'sim' set to TRUE when the length of the data is a multiple > 1 times the length of the positions")
			sim = TRUE
			}
		}
	else if(((lx/lcoords) %% sim)!=0){
		stop("The lengths of the data and positions do not agree")
		}
	# Get the logical version of 'sim', corresponsing to one set of positions for each time step. If sim=2 for example, the positions are used twice in ksmooth.SG():
	siml = sim==1
	
	# Get the original dimensions of the input:
	olddimcoords = dim(data$psxx)
	olddimx = dim(data[[var]])
	# Get the full dimensions of the acoustic data:
	D1 = max(data$lenb)
	D23 = data$numb[1]
	D3 = length(unique(c(data$freq)))
	D2 = D23 / D3
	D4 = length(data[[var]]) / (D1*data$numb[1])
	fulldim = c(D1, D2, D3, D4)
	#fulldim = c(max(data$lenb), data$numb[1]/length(unique(c(data$freq))), length(unique(data$freq)), length(data[[var]])/(max(data$lenb)*data$numb))
	pingdim = fulldim[seq_len(length(fulldim)-1)]
	timedim = tail(fulldim, 1)
	
	# Expand the lengths of 'h' and 'w':
	h = rep(h, length.out=length(pingdim))
	w = rep(w, length.out=length(pingdim))
	
	# Expand the dimensions to 3-D:
	# If 'sim' is FALSE, the position data has the same length as the acoustic data, and should also have the same dimensions:
	if(sim==0){
		dim(data$psxx) = fulldim
		dim(data$psyx) = fulldim
		dim(data$pszx) = fulldim
		}
	# If 'sim' is TRUE, the position data are only given for one time step, and should have the same dimension as the acoustic data at one single time step. The position data are recycled in ksmooth.SG() to be used in all time steps:
	else if(sim==1){
		dim(data$psxx) = pingdim
		dim(data$psyx) = pingdim
		dim(data$pszx) = pingdim
		}
	# If the position data should be rerecycled, but are given for more than one time step:
	else if(sim>0 && sim%%1==0){
		dim(data$psxx) = c(pingdim, timedim/sim)
		dim(data$psyx) = c(pingdim, timedim/sim)
		dim(data$pszx) = c(pingdim, timedim/sim)
		}
	dim(data[[var]]) = fulldim
		
	
	# If specified, drop redundant dimensions both in the data and in 'h' and 'w':
	if(drop){
		# Get the dimensions to drop:
		keepdim = which(pingdim!=1)
		# Drop elements from 'h' and 'w':
		h = h[keepdim]
		w = w[keepdim]
		
		if(sim==0){
			droppeddimpos = c(pingdim[keepdim], timedim)
			}
		else if(sim==1){
			droppeddimpos = pingdim[keepdim]
			}
		else if(sim>0 && sim%%1==0){
			droppeddimpos = c(pingdim[keepdim], timedim/sim)
			}
		droppeddimx = c(fulldim[keepdim], timedim)
		# Drop dimensions:
		dim(data$psxx) = droppeddimpos
		dim(data$psyx) = droppeddimpos
		dim(data$pszx) = droppeddimpos
		dim(data[[var]]) = droppeddimx
		}
	
	
	##### Execution #####
	data[[var]] = ksmooth.SG(data[c("psxx", "psyx", "pszx")], x=data[[var]], h=h, w=w, sim=sim, ind=ind, na.rm=na.rm)$x
	
	# Restore the original dimensions:
	dim(data$psxx) = olddimcoords
	dim(data$psyx) = olddimcoords
	dim(data$pszx) = olddimcoords
	dim(data[[var]]) = olddimx
	
	
	##### Output #####
	data
	##################################################
	##################################################
	}
