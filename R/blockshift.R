#*********************************************
#*********************************************
#' Returns a matrix with blocks of dimension 'd' stacked vertically or lined up horizontally depending on 'd' and the dimensions of 'x'. For a 3x16 matrix and d=c(3,4) the output will be a stack of 3x4 matrices constituting a 12x4 matrix. For a 16x3 matrix and d=c(4,3) the output will be a 4x12 matrix.
#'
#' @param x  is a matrix.
#' @param d  is the dimension of the blocks. If d[1]==nrow(x) blocks are shifted from horizontal to vertical structure, and if d[2]==ncol(x) blocks are shifted from vertical to horizontal structue.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname blockshift
#'
blockshift<-function(x,d){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-02-08 - Finished.
	########### DESCRIPTION: ###########
	# Returns a matrix with blocks of dimension 'd' stacked vertically or lined up horizontally depending on 'd' and the dimensions of 'x'. For a 3x16 matrix and d=c(3,4) the output will be a stack of 3x4 matrices constituting a 12x4 matrix. For a 16x3 matrix and d=c(4,3) the output will be a 4x12 matrix.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is a matrix.
	# ---d--- is the dimension of the blocks. If d[1]==nrow(x) blocks are shifted from horizontal to vertical structure, and if d[2]==ncol(x) blocks are shifted from vertical to horizontal structue.
	
	
	##################################################
	##################################################
	##### Preparation #####
	dimx=dim(x)
	
	
	##### Execution and output #####
	if(d[1]==dimx[1]){
		# Number of blocks:
		n=dimx[2]/d[2]
		if(n%%1 != 0){
			warning("Dimension mismatch, output is cropped")
			}
		n=floor(n)
		# Calculating indexes shifting the elements of 'x':
		a12=rep(1:d[1],n) + rep(seq(0,by=prod(d),l=n),each=d[1])
		a3=seq(0,by=d[1],l=d[2])
		a=outer(a12,a3,"+")
		# Clear some memory:
		rm(a12,a3)
		# Flatten 'a':
		dim(a)=NULL
		# Reorganize 'x':
		x=x[a]
		rm(a)
		# Add dimensions to 'x':
		dim(x)=c(prod(dimx)/d[2],d[2])
		# Return:
		x
		}
	else if(d[2]==dimx[2]){
		# Number of blocks:
		n=dimx[1]/d[1]
		if(n%%1 != 0){
			warning("Dimension mismatch, output is cropped")
			}
		n=floor(n)
		# Calculating indexes shifting the elements of 'x':
		a12=rep(1:d[1],d[2]) + rep(seq(0,by=dimx[1],l=d[2]),each=d[1])
		a3=seq(0,by=d[1],l=n)
		a=outer(a12,a3,"+")
		# Clear some memory:
		rm(a12,a3)
		# Flatten 'a':
		dim(a)=NULL
		# Reorganize 'x':
		x=x[a]
		rm(a)
		# Add dimensions to 'x':
		dim(x)=c(d[1],prod(dimx)/d[1])
		# Return:
		x
		}
	else{
		stop("Invalid dimensions 'd'")
		}
	##################################################
	##################################################
	}
