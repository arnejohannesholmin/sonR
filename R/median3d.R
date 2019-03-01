#*********************************************
#*********************************************
#' Median along the dimension given by 'along'. The same as apply(x,(1:3)[-along],median,na.rm=TRUE), but some hundred times faster.
#'
#' @param x  is a 3d array to be median smoothed along the 'along' dimension.
#' @param along Integer. The median is taken along this dimension.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD dim_all ones
#'
#' @export
#' @rdname median3d
#'
median3d=function(x, along=1){
	
	############### LOG: ###############
	# Start: 2012-01-24 - Clean version.

	d=dim_all(x)
	nd=length(d)
	add=3-nd
	if(nd<3){
		dim(x)=c(d,ones(add))
		d=dim_all(x)
		}
	else if(nd>3){
		stop("Only arrays up to three dimensions supported")
		}
	# Execution:
	U <- .C("median3d", 
		as.double(x), # (1) The data
		as.integer(d[1]), # (2) First dimension
		as.integer(d[2]), # (3) Second dimension
		as.integer(d[3]), # (4) Third dimension
		as.integer(along), # (5) The "along" dimension
		double(prod(d[-along])), # (6) The output data
		double(d[along]), # (7) The sorted values along the dimension given by along, used inside c++ to avoid memory limit
		NAOK=TRUE,
		PACKAGE="sonR")
	x = U[[6]]
	
	
	# Reset dimension to the original:
	if(add>0){
		dim(x)=NULL
		}
	else{
		dim(x)=d[-along]
		}
	
		
	x
}
