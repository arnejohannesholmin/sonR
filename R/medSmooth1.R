#*********************************************
#*********************************************
#' Median smoothing an array along the first dimension.
#'
#' @param x				An array to be median smoothed along the first dimension.
#' @param w				The width of median filter (odd number).
#' @param try.runmed	Locigal: If TRUE try the function \code{\link{medSmooth1_runmed}} first and then \code{\link{medSmooth1_median3d}}.
#' @param silent		Used in \code{\link{medSmooth1_runmed}}.
#' @param ...			Used in \code{\link{medSmooth1_runmed}} and \code{\link{medSmooth1_median3d}}.
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
	
	############### LOG: ###############
	# Start: 2013-09-19 - Clean version.
	
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
}
