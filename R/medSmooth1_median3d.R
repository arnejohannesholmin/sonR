#*********************************************
#*********************************************
#' Median smoothing an array along the first dimension.
#'
#' @param x  is an array to be median smoothed along the first dimension.
#' @param w  is the width of median filter (odd number).
#' @param maxsize  is the maximum size of the array used in the median smoothing with median3d(). If prod(dim(x)[1]+w-1,dim(x)[-1],w)>maxsize, an error is issued.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD NAs odd
#'
#' @export
#' @rdname medSmooth1_median3d
#'
medSmooth1_median3d<-function(x, w=3, maxsize=2e8){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-09-19 - Clean version.
	# Last: 2014-09-26 - Changed to use median3d of an expanded array due to bugs in runmed with NAs in the data.
	########### DESCRIPTION: ###########
	# Median smoothing an array along the first dimension.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is an array to be median smoothed along the first dimension.
	# ---w--- is the width of median filter (odd number).
	# ---maxsize--- is the maximum size of the array used in the median smoothing with median3d(). If prod(dim(x)[1]+w-1,dim(x)[-1],w)>maxsize, an error is issued.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	if(length(x)==0 || is.character(x)){
		return(x)
		}
	if(!odd(w)){
		warning(paste("'k' must be odd!  Changing 'k' to", w+1))
		w <- w+1
		}
	add <- (w-1)/2
	# Store the dimension of 'x', and convert to matrix:
	olddim <- dim(x)
	if(length(olddim)==0){
		dim(x) <- c(length(x),1)
		}
	else if(length(olddim)>2){
		dim(x) <- c(olddim[1],prod(olddim[-1]))
		}
	tempdim <- dim(x)
	tempdimLarge <- c(tempdim,w)
	tempdimLarge[1] <- tempdimLarge[1] + (w-1)
	
	
	########## Execution ##########
	# Create a 3D-array of dimension [w,tempdim], on which median3d() is run:
	if(prod(tempdimLarge)>maxsize){
		stop(paste0("An array of size ",prod(tempdimLarge)," elements will be created. To run, increase 'maxsize'"))
		}
	out <- rep(list(x),w)
	for(i in seq(0,w-1)){
		out[[i+1]] <- rbind(NAs(i,tempdim[2]), out[[i+1]], NAs(w-i-1,tempdim[2]))
		}
	out <- unlist(out, use.names=FALSE)
	dim(out) <- tempdimLarge
	# Or
	#out <- NAs(tempdimLarge)
	#for(i in seq_len(w)){
	#	out[seq(i,tempdim[1]+i-1),,i] <- x
	#	}
	out <- median3d(out[seq(add+1,tempdim[1]+add),,,drop=FALSE], along=3)
	
		
	########## Output ##########
	# Return the smoothed data with the original dimension:
	dim(out) <- olddim
	out
	##################################################
	##################################################
	}
