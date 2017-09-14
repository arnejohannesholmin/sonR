#*********************************************
#*********************************************
#' Applies a Gaussian 1-D, 2-D, or 3-D kernel on the data 'x' associated to locations 'coords', given the bandwidths given in the vector 'h'. When sim==1, NAs in the coordinates or the data are not supported.
#'
#' @param coords  is a list or matrix of spatial coordinates (with names "x", "y", and "z" for 3-D data, otherwise interpreted in that order).
#' @param x  is the data, possibly with one extra dimension in the case that sim==1, in which case the smoothing is done simultaneously for all steps in the last dimension of 'x' using the same coordinates 'coords'. If given as a vector (dim(x)=0), the length of 'x' must be an integer multiple of the length of the corrdinates, and if this integer is larger than 1, 'sim' is set to TRUE.
#' @param h  is the vector of bandwidths.
#' @param w  is the boundary of the kernel, outside which it is 0 (should have the same length as h).
#' @param sim  is a TRUE if smoothing should be done only along the first dimensions, simultaneously over the stages of the last dimension. If 'sim' is an integer larger than 1, the positions 'coords' are used 'sim' times, and the data 'x' should have length 'sim' times the length of one coordinate of 'coords'.
#' @param ind  is a list of subscripts for the dimensions of the data 'x', given in any form accepted by [], identifying where to smooth the data. Each element of the list corresponds to one dimension of 'x', and if 'ind' is shorter than the number of dimensions of 'x', no subsetting is done in the dimensions not given. 'ind' is treated by ind.expand().
#' @param na.rm  is single integer representing the dimension along which NAs are discarded from the smoothing in the case that sim==1. For example, if na.rm=2 and the dimension of 'x' is [5,12,7], and x[3,2:4,5]=NA, then all data x[,2:4,] will be excluded from the smoothing and set to NA. If na.rm=FALSE, no NAs should be contained in the data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD dim_all ind.expand
#' @importFrom utils tail
#'
#'
#' @export
#' @rdname ksmooth.SG
#'
ksmooth.SG=function(coords,x,h,w=h*3,sim=FALSE,ind=list(),na.rm=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-08-23 - Clean version.
	# Last: 2013-09-12 - Expanded to 1, 2, and 3 dimensions.
	########### DESCRIPTION: ###########
	# Applies a Gaussian 1-D, 2-D, or 3-D kernel on the data 'x' associated to locations 'coords', given the bandwidths given in the vector 'h'. When sim==1, NAs in the coordinates or the data are not supported.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---coords--- is a list or matrix of spatial coordinates (with names "x", "y", and "z" for 3-D data, otherwise interpreted in that order).
	# ---x--- is the data, possibly with one extra dimension in the case that sim==1, in which case the smoothing is done simultaneously for all steps in the last dimension of 'x' using the same coordinates 'coords'. If given as a vector (dim(x)=0), the length of 'x' must be an integer multiple of the length of the corrdinates, and if this integer is larger than 1, 'sim' is set to TRUE.
	# ---h--- is the vector of bandwidths.
	# ---w--- is the boundary of the kernel, outside which it is 0 (should have the same length as h).
	# ---sim--- is a TRUE if smoothing should be done only along the first dimensions, simultaneously over the stages of the last dimension. If 'sim' is an integer larger than 1, the positions 'coords' are used 'sim' times, and the data 'x' should have length 'sim' times the length of one coordinate of 'coords'.
	# ---ind--- is a list of subscripts for the dimensions of the data 'x', given in any form accepted by [], identifying where to smooth the data. Each element of the list corresponds to one dimension of 'x', and if 'ind' is shorter than the number of dimensions of 'x', no subsetting is done in the dimensions not given. 'ind' is treated by ind.expand().
	# ---na.rm--- is single integer representing the dimension along which NAs are discarded from the smoothing in the case that sim==1. For example, if na.rm=2 and the dimension of 'x' is [5,12,7], and x[3,2:4,5]=NA, then all data x[,2:4,] will be excluded from the smoothing and set to NA. If na.rm=FALSE, no NAs should be contained in the data.
	

	##################################################
	##################################################
	##### Preparation #####
	# Get the dimension of the data 'x':
	DIM=length(dim(x))
	# If 'x' is given as a vector, try to expand it to an array with dimensions fitting the dimensions of the 'coords':
	if(DIM==0){
		if(is.list(coords)){
			numt=length(x)/prod(sapply(coords,length))
			DIM=length(coords)
			if(numt==1){
				dim(x)=c(length(x),rep(1,DIM-1))
				}
			else if(numt%%1==0){
				dim(x)=c(length(x),rep(1,DIM-1,numt))
				sim=TRUE
				}
			else{
				stop("The dimensions of 'x' do not fit the dimensions of the list 'coords'")
				}
			}
		else{
			numt=length(x)/NROW(coords)
			DIM=NCOL(coords)
			if(numt==1){
				dim(x)=c(length(x),rep(1,DIM-1))
				}
			else if(numt%%1==0){
				dim(x)=c(length(x),rep(1,DIM-1,numt))
				sim=TRUE
				}
			else{
				stop("The dimensions of 'x' do not fit the dimensions of the matrix 'coords'")
				}
			}
		}
	siml=sim==1
	# 'sim' is here interpreted as the number of repitisions of the position data, and if sim==FALSE, 'sim' is reset to 1:
	if(sim==0){
		sim=1
		}
		
	# If sim==1, smoothing should not be performed over the last dimension of 'x':
	DIM=DIM-siml
	
	# Define the valid coordinate names:
	validCoords=c("x","y","z")[seq_len(DIM)]
	# Get the dimension and order of the coordinates 'coords':
	if(is.list(coords)){
		# Identify the elements of 'coords' which has names "x", "y" or "z":
		withNames=names(coords) %in% validCoords
		# Add any wissing coordinate names:
		if(any(!withNames)){
			missingCoordsNames=setdiff(validCoords,names(coords)[withNames])
			nTooMany=sum(!withNames)-length(missingCoordsNames)
			if(nTooMany>0){
				names(coords)[!withNames]=c(missingCoordsNames,paste(rep("z",nTooMany),seq_len(nTooMany),sep=""))
				}
			}
		# Get the order of the coordinates:
		orderCoords=order(names(coords))
		ndimCoords=min(4,length(dim(coords[[1]])))
		lengthCoords=lapply(coords,length)
		}
	else{
		if(length(coords)>0){
			dimcoords=dim(coords)
			# If only one vector is given:
			if(length(dimcoords)==0){
				dim(coords)=c(length(coords),1)
				orderCoords=1
				ndimCoords=1
				}
			else{
				# Identify the elements of 'coords' which has names "x", "y" or "z":
				withNames=colnames(coords) %in% validCoords
				# Add any wissing coordinate names:
				if(any(!withNames)){
					missingCoordsNames=setdiff(validCoords,colnames(coords)[withNames])
					nTooMany=sum(!withNames)-length(missingCoordsNames)
					if(nTooMany>0){
						colnames(coords)[!withNames]=c(missingCoordsNames,paste(rep("z",nTooMany),seq_len(nTooMany),sep=""))
						}
					}	
				# Get the order of the coordinates:
				orderCoords=order(colnames(coords))
				ndimCoords=min(4,ncol(coords))
				}
			# Expand the dimensions of 'coords' to those of the data 'x':
			if(prod(dim(x)[seq_len(DIM)])/nrow(coords) == sim){
				dim(coords)=c(dim(x)[seq_len(DIM-1)],dim(x)[DIM]/sim,ndimCoords)
				}	
			else{
				stop("The length of the coordinate vectors 'coords' must equal the length of the data 'x'")
				}
			}
		}
	
	# Test for agreement between the 'coords' and the 'x':
	if(DIM!=ndimCoords){
		stop("The number of coordinates and the number of dimensions of the data 'x' do not agree. If sim==1, 'x' is allowed to have one dimension higher than the number of coordinates. Otherwise, the number of dimensions of 'x' and the number of coordinates should be equal")
		}
	# Test for too low dimensions:
	if(DIM==0){
		warning("When the data 'x' has only one dimension, sim==1 is not valid")
		return(list(coords=coords, x=x))
		}
	# Test for high dimensions:
	if(DIM==4){
		warning("Four-dimensional data are smoothed one time step at the time in a for loop")
		}
	else if(DIM>4){
		warning("Only 1-D, 2-D, or 3-D smoothing supported. For smoothing 4-D data along the first three dimensions simultaneously over the stages of the fourth dimension, use sim==1")
		return(list(coords=coords, x=x))
		}
	dimx=dim_all(x)
	
	# Check the lengths of 'h' and 'w':
	if(length(h)<DIM){
		h=rep(h,length.out=DIM)
		}
	if(length(w)<DIM){
		w=rep(w,length.out=DIM)
		}
	
	# If all elements of 'h' are equal to 0, retunr the unsmoothed data:
	if(all(h==0)){
		return(list(coords=coords, x=x, h=h, w=w, sim=sim))
		}
	
	# Treat the subset-indices defining where to smooth:
	ind=ind.expand(ind,dimx)
	
	
	##### Execution #####
	# Apply the spatial kernel smoother with simultaneous smoothing for all stages of the last dimension:
	if(siml){
		
		# If the data has 2 dimensions and the coordinates have 1:
		if(DIM==1){
			cat("Smoothing in 1D simultaneously over all time steps\n")
			# Detect NAs and discard from the smoothing:
			if(na.rm==1){
				nas=which(rowSums(is.na(x))>0)
				x[nas,]=NA
				ind[[1]]=setdiff(ind[[1]],nas)
				}
			# Coordinates in a list:
			if(is.list(coords)){
				x[ ind[[1]], ind[[2]] ] <- .C( "kernSmooth1dGaussMultipleNoNA", 
				# The positions:
				as.double(coords[[orderCoords[1]]][ ind[[1]] ]), 
				# The data:
				as.double(x[ ind[[1]], ind[[2]] ]), 
				# The dimensions of the data:
				as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), 
				# The smoothing bandwidth:
				as.double(h[1]), 
				# The boundary of the kernel:
				as.double(w[1]), 
				# The output data:
				as.double(x[ ind[[1]], ind[[2]] ]), PACKAGE="sonR")[[7]]
				}	
			# Coordinates in a matrix:
			else{
				x[ ind[[1]], ind[[2]] ] <- .C( "kernSmooth1dGaussMultipleNoNA", 
				# The positions:
				as.double(coords[ind[[1]],orderCoords[1]]), 
				# The data:
				as.double(x[ ind[[1]], ind[[2]] ]), 
				# The dimensions of the data:
				as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), 
				# The smoothing bandwidth:
				as.double(h[1]), 
				# The boundary of the kernel:
				as.double(w[1]), 
				# The output data:
				as.double(x[ ind[[1]], ind[[2]] ]), PACKAGE="sonR")[[7]]
				}
			}
		
		# If the data has 3 dimensions and the coordinates have 2:
		else if(DIM==2){
			cat("Smoothing in 2D simultaneously over all time steps\n")
			# Detect NAs and discard from the smoothing:
			if(na.rm==1){
				nas=which(rowSums(is.na(x))>0)
				x[nas,,]=NA
				ind[[1]]=setdiff(ind[[1]],nas)
				}
			else if(na.rm==2){
				nas=which(rowSums(colSums(is.na(x)))>0)
				x[,nas,]=NA
				ind[[2]]=setdiff(ind[[2]],nas)
				}
			# Coordinates in a list:
			if(is.list(coords)){
				x[ ind[[1]], ind[[2]], ind[[3]] ] <- .C( "kernSmooth2dGaussMultipleNoNA", 
				# The positions:
				as.double(coords[[orderCoords[1]]][ ind[[1]], ind[[2]] ]), as.double(coords[[orderCoords[2]]][ ind[[1]], ind[[2]] ]), 
				# The data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]] ]), 
				# The dimensions of the data:
				as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(ind[[3]])), 
				# The smoothing bandwidth:
				as.double(h[1]), as.double(h[2]), 
				# The boundary of the kernel:
				as.double(w[1]), as.double(w[2]), 
				# The output data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]] ]), PACKAGE="sonR")[[11]]
				}	
			# Coordinates in a matrix:
			else{
				x[ ind[[1]], ind[[2]], ind[[3]] ] <- .C( "kernSmooth2dGaussMultipleNoNA", 
				# The positions:
				as.double(coords[ind[[1]],ind[[2]],orderCoords[1]]), as.double(coords[ind[[1]],ind[[2]],orderCoords[2]]), 
				# The data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]] ]), 
				# The dimensions of the data:
				as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(ind[[3]])), 
				# The smoothing bandwidth:
				as.double(h[1]), as.double(h[2]), 
				# The boundary of the kernel:
				as.double(w[1]), as.double(w[2]), 
				# The output data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]] ]), PACKAGE="sonR")[[11]]
				}
			}
		
		# If the data has 4 dimensions and the coordinates have 3:
		else if(DIM==3){
			cat("Smoothing in 3D simultaneously over all time steps\n")
			# Detect NAs and discard from the smoothing:
			if(na.rm==1){
				nas=which(rowSums(is.na(x))>0)
				x[nas,,,]=NA
				ind[[1]]=setdiff(ind[[1]],nas)
				}
			else if(na.rm==2){
				nas=which(rowSums(colSums(is.na(x)))>0)
				x[,nas,,]=NA
				ind[[2]]=setdiff(ind[[2]],nas)
				}
			else if(na.rm==3){
				nas=which(rowSums(colSums(is.na(x)))>0)
				x[,,nas,]=NA
				ind[[3]]=setdiff(ind[[3]],nas)
				}
			# Coordinates in a list:
			if(is.list(coords)){
				x[ ind[[1]], ind[[2]], ind[[3]], ind[[4]] ] <- .C( "kernSmooth3dGaussMultipleNoNA", 
				# The positions:
				as.double(coords[[orderCoords[1]]][ ind[[1]], ind[[2]], ind[[3]] ]), as.double(coords[[orderCoords[2]]][ ind[[1]], ind[[2]], ind[[3]] ]), as.double(coords[[orderCoords[3]]][ ind[[1]], ind[[2]], ind[[3]] ]), 
				# The data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]], ind[[4]] ]), 
				# The dimensions of the data:
				as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(ind[[3]])), as.integer(length(ind[[4]])), 
				# The smoothing bandwidth:
				as.double(h[1]), as.double(h[2]), as.double(h[3]), 
				# The boundary of the kernel:
				as.double(w[1]), as.double(w[2]), as.double(w[3]), 
				# The output data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]], ind[[4]] ]), PACKAGE="sonR")[[15]]
				}	
			# Coordinates in a matrix:
			else{
				x[ ind[[1]], ind[[2]], ind[[3]], ind[[4]] ] <- .C( "kernSmooth3dGaussMultipleNoNA", 
				# The positions:
				as.double(coords[ind[[1]],ind[[2]],ind[[3]],orderCoords[1]]), as.double(coords[ind[[1]],ind[[2]],ind[[3]],orderCoords[2]]), as.double(coords[ind[[1]],ind[[2]],ind[[3]],orderCoords[3]]), 
				# The data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]], ind[[4]] ]), 
				# The dimensions of the data:
				as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(ind[[3]])), as.integer(length(ind[[4]])), 
				# The smoothing bandwidth:
				as.double(h[1]), as.double(h[2]), as.double(h[3]), 
				# The boundary of the kernel:
				as.double(w[1]), as.double(w[2]), as.double(w[3]), 
				# The output data:
				as.double(x[ ind[[1]], ind[[2]], ind[[3]], ind[[4]] ]) , PACKAGE="sonR")[[15]]
				}
			}
		}
	
	# Apply the standard spatial kernel smoother, allowing for missing values in the data:
	else{
		x[is.na(x)] = NaN
		
		# Get the number of time steps as the last dimension of 'x':
		nt=tail(dimx,1)/sim
		if(nt %% 1 != 0){
			stop("When 'sim' is an integer number of times to recycle the position data 'coord', the data 'x' should have a number of time steps equal to a multiple 'sim' times the number of time steps of the positions")
			}

		for(i in seq_len(sim)){
			# Get the time step indices for step 'i' in the for loop:
			timei=1:nt + (i-1)*nt
			
			# If the data and coordinates have 1 dimension:
			if(DIM==1){
				# Get the time steps present in the last element of 'ind', and select the subset 'valid' for the positions and 'timei[valid]' for the data:
				valid=which(timei %in% ind[[1]])
				cat("Smoothing in 1D\n")
				
				if(length(valid)>0){
					# Coordinates in a list:
					if(is.list(coords)){
						x[ timei[valid] ] <- .C( "kernSmooth1dGauss", 
						# The positions:
						as.double(coords[[orderCoords[1]]][ valid ]), 
						# The data:
						as.double(x[ timei[valid] ]), 
						# The dimensions of the data:
						as.integer(length(valid)), 
						# The smoothing bandwidth:
						as.double(h[1]), 
						# The boundary of the kernel:
						as.double(w[1]), 
						# The output data:
						as.double(x[ timei[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[6]]
						}	
					# Coordinates in a matrix:
					else{
						x[ timei[valid] ] <- .C( "kernSmooth1dGauss", 
						# The positions:
						as.double(coords[ valid, orderCoords[1] ]), 
						# The data:
						as.double(x[ timei[valid] ]), 
						# The dimensions of the data:
						as.integer(length(valid)), 
						# The smoothing bandwidth:
						as.double(h[1]),  
						# The boundary of the kernel:
						as.double(w[1]), 
						# The output data:
						as.double(x[ timei[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[6]]
						}
					}
				}
			
			# If the data and coordinates have 2 dimensions:
			else if(DIM==2){
				# Get the time steps present in the last element of 'ind', and select the subset 'valid' for the positions and 'timei[valid]' for the data:
				valid=which(timei %in% ind[[2]])
				cat("Smoothing in 2D\n")
				
				if(length(valid)>0){
					# Coordinates in a list:
					if(is.list(coords)){
						x[ ind[[1]], timei[valid] ] <- .C( "kernSmooth2dGauss", 
						# The positions:
						as.double(coords[[orderCoords[1]]][ ind[[1]], valid ]), as.double(coords[[orderCoords[2]]][ ind[[1]], valid ]), 
						# The data:
						as.double(x[ ind[[1]], timei[valid] ]), 
						# The dimensions of the data:
						as.integer(length(ind[[1]])), as.integer(length(valid)), 
						# The smoothing bandwidth:
						as.double(h[1]), as.double(h[2]), 
						# The boundary of the kernel:
						as.double(w[1]), as.double(w[2]), 
						# The output data:
						as.double(x[ ind[[1]], timei[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[10]]
						}	
					# Coordinates in a matrix:
					else{
						x[ ind[[1]], timei[valid] ] <- .C( "kernSmooth2dGauss", 
						# The positions:
						as.double(coords[ ind[[1]], valid, orderCoords[1] ]), as.double(coords[ ind[[1]], valid, orderCoords[2] ]), 
						# The data:
						as.double(x[ ind[[1]], timei[valid] ]), 
						# The dimensions of the data:
						as.integer(length(ind[[1]])), as.integer(length(valid)), 
						# The smoothing bandwidth:
						as.double(h[1]), as.double(h[2]), 
						# The boundary of the kernel:
						as.double(w[1]), as.double(w[2]), 
						# The output data:
						as.double(x[ ind[[1]], timei[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[10]]
						}
					}
				}
			
			# If the data and coordinates have 3 dimension:
			else if(DIM==3){
				# Get the time steps present in the last element of 'ind', and select the subset 'valid' for the positions and 'timei[valid]' for the data:
				valid=which(timei %in% ind[[3]])
				cat("Smoothing in 3D\n")
				
				if(length(valid)>0){
					# Coordinates in a list:
					if(is.list(coords)){
						x[ ind[[1]], ind[[2]], timei[valid] ] <- .C( "kernSmooth3dGauss", 
						# The positions:
						as.double(coords[[orderCoords[1]]][ ind[[1]], ind[[2]], valid ]), as.double(coords[[orderCoords[2]]][ ind[[1]], ind[[2]], valid ]), as.double(coords[[orderCoords[3]]][ ind[[1]], ind[[2]], valid ]), 
						# The data:
						as.double(x[ ind[[1]], ind[[2]], timei[valid] ]), 
						# The dimensions of the data:
						as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(valid)), 
						# The smoothing bandwidth:
						as.double(h[1]), as.double(h[2]), as.double(h[3]), 
						# The boundary of the kernel:
						as.double(w[1]), as.double(w[2]), as.double(w[3]), 
						# The output data:
						as.double(x[ ind[[1]], ind[[2]], timei[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[14]]
						}	
					# Coordinates in a matrix:
					else{
						x[ ind[[1]], ind[[2]], timei[valid] ] <- .C( "kernSmooth3dGauss", 
						# The positions:
						as.double(coords[ ind[[1]], ind[[2]], valid, orderCoords[1] ]), as.double(coords[ ind[[1]], ind[[2]], valid, orderCoords[2] ]), as.double(coords[ ind[[1]], ind[[2]], valid, orderCoords[3] ]), 
						# The data:
						as.double(x[ ind[[1]], ind[[2]], timei[valid] ]), 
						# The dimensions of the data:
						as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(valid)), 
						# The smoothing bandwidth:
						as.double(h[1]), as.double(h[2]), as.double(h[3]), 
						# The boundary of the kernel:
						as.double(w[1]), as.double(w[2]), as.double(w[3]), 
						# The output data:
						as.double(x[ ind[[1]], ind[[2]], timei[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[14]]
						}
					}
				}
			else if(DIM==4){
				cat("Smoothing in 3D for each time step\n")
				# Move through the time steps of the current  step 'i' in the for loop:
				for(j in seq_len(nt)){
					# Get the time step indices for step 'i' in the for loop:
					timej=j + (i-1)*nt
			
					# Coordinates in a list:
					if(is.list(coords)){
						# Get the time steps present in the last element of 'ind', and select the subset 'valid' for the positions and 'timei[valid]' for the data:
						valid=which(timej %in% ind[[4]])
						
						if(length(valid)>0){
							x[ ind[[1]], ind[[2]], ind[[3]], timej[valid] ] <- .C( "kernSmooth3dGauss", 
							# The positions:
							as.double(coords[[orderCoords[1]]][ ind[[1]], ind[[2]], ind[[3]], valid ]), as.double(coords[[orderCoords[2]]][ ind[[1]], ind[[2]], ind[[3]], valid ]), as.double(coords[[orderCoords[3]]][ ind[[1]], ind[[2]], ind[[3]], valid ]), 
							# The data:
							as.double(x[ ind[[1]], ind[[2]], ind[[3]], timej[valid] ]), 
							# The dimensions of the data:
							as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(ind[[3]])), 
							# The smoothing bandwidth:
							as.double(h[1]), as.double(h[2]), as.double(h[3]), 
							# The boundary of the kernel:
							as.double(w[1]), as.double(w[2]), as.double(w[3]), 
							# The output data:
							as.double(x[ ind[[1]], ind[[2]], ind[[3]], timej[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[14]]
							}	
						# Coordinates in a matrix:
						else{
							x[ ind[[1]], ind[[2]], ind[[3]], timej[valid] ] <- .C( "kernSmooth3dGauss", 
							# The positions:
							as.double(coords[ ind[[1]], ind[[2]], ind[[3]], valid, orderCoords[1] ]), as.double(coords[ ind[[1]], ind[[2]], ind[[3]], valid, orderCoords[2] ]), as.double(coords[ ind[[1]], ind[[2]], ind[[3]], valid, orderCoords[3] ]), 
							# The data:
							as.double(x[ ind[[1]], ind[[2]], ind[[3]], timej[valid] ]), 
							# The dimensions of the data:
							as.integer(length(ind[[1]])), as.integer(length(ind[[2]])), as.integer(length(ind[[3]])), 
							# The smoothing bandwidth:
							as.double(h[1]), as.double(h[2]), as.double(h[3]), 
							# The boundary of the kernel:
							as.double(w[1]), as.double(w[2]), as.double(w[3]), 
							# The output data:
							as.double(x[ ind[[1]], ind[[2]], ind[[3]], timej[valid] ]), NAOK=TRUE, PACKAGE="sonR")[[14]]
							}
						}
					}
				}
			}
		x[is.na(x)] = NA
		}
	
	# Reset the dimension of 'x':
	dim(x)=dimx
	# Restore dimensions of 'corrds' if initially given as a matrix:
	if(!is.list(coords)){
		dim(coords)=c(length(coords)/DIM,DIM)
		}
	
	
	##### Output #####
	# Return a list of the coordinates, data, and parameters:
	list(coords=coords, x=x, h=h, w=w, sim=sim)
	##################################################
	##################################################
	}
