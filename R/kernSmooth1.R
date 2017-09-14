#*********************************************
#*********************************************
#' (Usually Gaussian) smoothing an array along the first dimension. Note: Due to high variability in the speed of the fft() used in convolve() in this function, the performance times vary greatly for different lengths of the input data. NAs are removed in the actual convolution, so a large number of NAs should not slow down the function.
#'
#' @param x  is an array to be smoothed along the first dimension.
#' @param kern  is the kernel of the smoothing filter, defaulted to a Gaussian kernel.
#' @param nsd  is the number of standard deviations on either side of the mean of the Gaussian kernel, outside which the kernel is zero.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom stats convolve
#'
#' @export
#' @rdname kernSmooth1
#'
kernSmooth1<-function(x, kern=3, nsd=3){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-11-23 - Clean version.
	########### DESCRIPTION: ###########
	# (Usually Gaussian) smoothing an array along the first dimension.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is an array to be smoothed along the first dimension.
	# ---kern--- is the kernel of the smoothing filter, defaulted to a Gaussian kernel.
	# ---nsd--- is the number of standard deviations on either side of the mean of the Gaussian kernel, outside which the kernel is zero.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# If the kernel is given as a single number, it is interpreted as the standard deviation of a Gaussian kernel:
	if(length(kern)==1){
		kern <- dnorm(seq(floor(-kern*nsd), ceiling(kern*nsd)), sd=kern)
		kern <- kern/sum(kern)
		}
	# Store the dimension of 'x', and convert to matrix:
	olddim <- dim(x)
	if(length(olddim)==0){
		dim(x) <- c(length(x),1)
		}
	else if(length(olddim)>2){
		dim(x) <- c(olddim[1], prod(olddim[-1]))
		}
	numb <- dim(x)[2]
	lenb <- dim(x)[1]
	keep <- (length(kern)+1)/2
	
	
	########## Execution ##########
	# Run through the columns and smooth:
	for(i in seq_len(numb)){
		
		# Try to replace NAs by 0 and convolve first the array, and then the array where all positive values are replaced by 1, and divide the first result by the second result:
		atNA <- is.na(x[,i])
		x0 <- replace(x[,i], atNA, 0)
		x01 <- replace(x0, !atNA, 1)
		smooth <- convolve(x0, kern, type="o")[seq(keep,length.out=lenb)]
		scale <- convolve(x01, kern, type="o")[seq(keep,length.out=lenb)]
		x[,i] <- replace(smooth / scale, atNA, NA)
		# 
		# 
		# 
		# whichnotNA <- which(!is.na(x[,i]))
		# if(any(diff(whichnotNA)>1)){
		# 	warning("Missing values should only be located at the beginning or end of the vectors")
		# }
		# if(length(whichnotNA)){
		# 	# Apply the kernel:
		# 	x[whichnotNA,i] <- convolve(x[whichnotNA,i], kern, type="o")[seq(keep,length.out=length(whichnotNA))]
		# 	# Scale up where the entire kernel has not been convolved:
		# 	lkern <- length(kern)
		# 	lx <- length(whichnotNA)
		# 	s <- seq_len(lx)
		# 	m <- (lkern + 1) / 2
		# 	first <- pmax(1, 1 + m - s)
		# 	last <- pmin(lkern, rev(1 + lkern - (1 + m - s)))
		# 	firstlast <- cbind(first, last)
		# 	combinations <- apply(firstlast, 1, paste, collapse="-")
		# 	# Detect non-duplicated combinations, and get the sum of the kernel for these combinations, and divide by these values:
		# 	notdup <- which(!duplicated(combinations))
		# 	kernsums <- sapply(notdup, function(i) sum(kern[seq(first[i], last[i])]))
		# 	names(kernsums) <- combinations[notdup]
		# 	x[whichnotNA,i] <- unname(x[whichnotNA,i] / kernsums[combinations])
		# }
		# At both ends, comvolve manually to avoid the effect of padding with zeros done by convolve():
		
		}
	
		
	########## Output ##########
	# Return the smoothed data with the original dimension:
	dim(x) <- olddim
	return(x)
	##################################################
	##################################################
	}
