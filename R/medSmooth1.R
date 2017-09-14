#*********************************************
#*********************************************
#' Median smoothing an array along the first dimension.
#'
#' @param x  is an array to be median smoothed along the first dimension.
#' @param w  is the width of median filter (odd number).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname medSmooth1
#'
medSmooth1<-function(x, w=3, try.runmed=FALSE, silent=FALSE, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-09-19 - Clean version.
	########### DESCRIPTION: ###########
	# Median smoothing an array along the first dimension.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is an array to be median smoothed along the first dimension.
	# ---w--- is the width of median filter (odd number).
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Try the fast method:
	m <- NULL
	if(try.runmed){
		m <- try(medSmooth1_runmed(x, w, ...), silent=silent)
	}
	# ... and if failing due to NAs, try the slow method:
	if(length(m)==0 || !is.numeric(m)){
		if(try.runmed){
			warning("medSmooth1_runmed() failed due to NAs")
		}
		m <- medSmooth1_median3d(x, w, ...)
	}
	return(m)
	##################################################
	##################################################
}
