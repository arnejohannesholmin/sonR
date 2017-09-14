#*********************************************
#*********************************************
#' Expands the fucntion filter() to the edges, instead of returning NA.
#'
#' @param x  is the data vector.
#' @param filter  is a vector of the filter values as input to filter().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD odd ones zeros
#' @importFrom stats filter
#'
#' @export
#' @rdname filter2
#'
filter2<-function(x,filter){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-06-13 - Clean version.
	########### DESCRIPTION: ###########
	# Expands the fucntion filter() to the edges, instead of returning NA.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is the data vector.
	# ---filter--- is a vector of the filter values as input to filter().
	
	
	##################################################
	##################################################
	# Different action for filters of odd and even length:
	lf=length(filter)
	# If odd, simply add (lf-1)/2 zeros on both sides of x, and remove the missing values in the output:
	if(odd(lf)){
		lz=(lf-1)/2
		# Pad the vector with zeros, and devide by the filtered flat response:
		z=zeros(lz)
		out=filter(c(z,x,z),filter)/filter(c(z,ones(length(x)),z),filter)
		out[lz+seq_along(x)]
		}
	else{
		lz=(lf)/2
		# Pad the vector with zeros, and devide by the filtered flat response:
		z1=zeros(lz-1)
		z2=zeros(lz)
		out=filter(c(z1,x,z2),filter)/filter(c(z1,ones(length(x)),z2),filter)
		out[lz-1+seq_along(x)]
		}
	##################################################
	##################################################
	}
