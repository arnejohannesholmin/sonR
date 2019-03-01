#*********************************************
#*********************************************
#' Median smoothing an array along the first dimension.
#'
#' @param x  is an array to be median smoothed along the first dimension.
#' @param w  is the width of median filter (odd number).
#' @param ...  Used for robustness.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD odd zeros
#' @importFrom stats runmed
#'
#' @export
#' @rdname medSmooth1_runmed
#'
medSmooth1_runmed<-function(x, w=3, ...){
	
	############### LOG: ###############
	# Start: 2013-09-19 - Clean version.
	
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
		dim(x) <- c(olddim[1], prod(olddim[-1]))
		}
	tempdim <- dim(x)
	
	
	########## Execution ##########
	# In order to handle NAs, add zeros at the end of the columns, collapse to a vector, apply runmed, restore the dimensions, and remove the zeros:
	x <- rbind(zeros(add,tempdim[2]), x, zeros(add,tempdim[2]))
	largedim <- dim(x)
	#x <- runmed(c(x),w,...)
	x <- runmed(c(x), w)
	dim(x) <- largedim
	x <- x[seq_len(tempdim[1]) + add,]
	# Run through the columns and smooth:
	#for(i in seq_len(ncol(x))){
	#	x[,i] <- runmed(x[,i],w,...)
	#	}
	
		
	########## Output ##########
	# Return the smoothed data with the original dimension:
	dim(x) <- olddim
	x
}
