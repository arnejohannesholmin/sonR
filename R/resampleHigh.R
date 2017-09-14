#*********************************************
#*********************************************
#' Resamples large values from the remaining values (useful in analyses containing suspucously high values which desrupt the analysis).
#'
#' @param x  is the vector to modify.
#' @param high  is the boundary classifying values as high values, defined by x/mean(x,trim=trim) (values larger than 'high' times the trimmed mean are defined as high values).
#' @param trim  is used in x/mean(x,trim=trim).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname resampleHigh
#'
resampleHigh=function(x,high,trim=0.1){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-12-05 - Clean version.
	########### DESCRIPTION: ###########
	# Resamples large values from the remaining values (useful in analyses containing suspucously high values which desrupt the analysis).
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is the vector to modify.
	# ---high--- is the boundary classifying values as high values, defined by x/mean(x,trim=trim) (values larger than 'high' times the trimmed mean are defined as high values).
	# ---trim--- is used in x/mean(x,trim=trim).
	
	
	##################################################
	##################################################
	toohigh=x/mean(x,trim=trim) > high
	if(sum(toohigh)>0){
		x[toohigh]=sample(x[!toohigh],sum(toohigh))
		}
	x
	##################################################
	##################################################
	}
