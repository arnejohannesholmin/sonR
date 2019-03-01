#*********************************************
#*********************************************
#' Extrapolates or shrinks a matrix to according to the specifications given by 'a', 'b' and 'along'. If 'along' (defaulted to 0) is not in 1:nrow(x), 'a' and 'b' are interpreted as the number of rows to add at the beginning and at the end of the columns of the matrix (negative values imply shrinkage). If along is in in 1:nrow(x), 'a' and 'b' are interpreted as the desired lower and upper value to which the values of x[,along] are to be extrapolated by a simple linear expansion of the two first and the two last rows of 'x'. The column given by 'along' should be sorted.
#'
#' @param x  is the input vector or matrix/data.frame. If 'x' is a vector, it is interpreted as a matrix/data.frame of one column
#' @param a,b give the range to which the column given by 'along' is extrapolated/shrinked. If along is not in 1:nrow(x), 'a' and 'b' are interpreted as the number of rows to add at the beginning and at the end of the columns (negative values imply shrinkage). If b==NULL and 'a' has length > 1, the second element of 'a' is interpreted as 'b'.
#' @param along  is the column along which the extrapolation is specified.
#' @param grid  is TRUE if the first and the last row of the output are on the grid defined by the expansion.
#' @param byrow  is TRUE if the expansion/shrinkage is to be done on the rows instead of the columns.
#' @param drop.out  is TRUE if the output should be stripped of empty dimensions by the drop() function.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ones
#'
#' @export
#' @rdname extrapolate.matrix
#'
extrapolate.matrix<-function(x, a=NULL, b=NULL, along=0, grid=FALSE, byrow=FALSE, drop.out=TRUE){
	
	############### LOG: ###############
	# Start: 2009-07-08 - Clean version.
		
	# If 'along' is a string and matches one of the column names of 'x', this column is chosen:
	if(is.character(along)){
		matchalong=match(along,colnames(x))
		if(!is.na(matchalong)){
			along=matchalong
			}
		else{
			warning("'along' does not match any of the colnames of 'x'. Default chosen")
			along=0
			}
		}
	
	if(is.data.frame(x)){
		x=as.matrix(x)
		}
		
	# Number of dimensions of 'x':
	ndimx=length(dim(x))
	# 'x' must be a matrix/data.frame or a vector of length >1:
	if(ndimx>2){
		stop("Higher dimensional inputs are not supported")
		}
	if(ndimx<2){
		lx=length(x)
		if(lx<2){
			return(rep(x,sum(a)))
			}
		else{
			dim(x)=c(lx,1)
			}
		}
		
	# If byrow==TRUE the rows are treated as columns:
	if(byrow){
		x=t(x)
		}
	
	# Remove duplicate rows at the start and beginning of the matrix:
	while(TRUE){
		if(all(x[1,]==x[2,])){
			x=x[-1,]
			}
		else{
			break
			}
		}
	while(TRUE){
		if(all(x[nrow(x),]==x[nrow(x)-1,])){
			x=x[-nrow(x),]
			}
		else{
			break
			}
		}
	
	# Dimensions of 'x':
	nrowx=nrow(x)
	ncolx=ncol(x)
	
	# If 'x' has only one row, interpolation cannot be done:
	if(nrowx==1){
		return(drop(x))
		}
	
	# 'a' and 'b' may both be given in 'a': 
	if(is.null(b) && length(a)>1){
		b=a[2]
		a=a[1]
		}
	else{
		a=a[1]
		b=b[1]
		}
	# If grid==FALSE and simple==FALSE, 'olda' and 'oldb' are used at the end of the function when adding the remainder of the matrix so that it reaches the range defiened by c(a,b):
	olda=a
	oldb=b
	
	# 'simple' is TRUE if 'a' and 'b' simply specifies the number of values to add at the beginning and the end of the columns:
	simple=!along %in% (1:ncolx)
	
	
	##### Execution #####
	if(!is.null(b)){
		# Defining the end row and the difference between the end row and the previous row:
		xn=x[nrowx,]
		diffxn=xn-x[nrowx-1,]
		# If 'a' and 'b' specifies the number of values to add at the beginning and the end of the columns:
		if(simple){
			b=floor(b)
			# Shrink:
			if(b<1){
				ind=(1:nrowx)<=(nrowx+b)
				# 'x' shrinked to nothing:
				if(!any(ind)){
					return(NULL)
					}
				else{
					x=x[ind,,drop=FALSE]
					}
				}
			# Extrapolate:
			else{
				x=rbind(x,outer(sequence(b),diffxn)+outer(ones(b),xn))
				}
			}
		# If 'a' and 'b' specifies the range to which the column given by 'along' is extrapolated/shrinked:
		else{
			toright=(b-xn[along])/diffxn[along]
			# Shrink:
			if(toright<1){
				ind=x[,along]<=b
				# 'x' shrinked to nothing:
				if(!any(ind)){
					return(NULL)
					}
				else{
					x=x[ind,,drop=FALSE]
					}
				}
			# Extrapolate:
			else{
				x=rbind(x,outer(sequence(toright),diffxn)+outer(ones(toright),xn))
				}
			}
		nrowx=nrow(x)
		}
		
	if(!is.null(a)){
		# Defining the end row and the difference between the end row and the previons row:
		x1=x[1,]
		diffx1=x1-x[2,]
		# If 'a' and 'b' specifies the number of values to add at the beginning and the end of the columns:
		if(simple){
			a=floor(a)
			# Shrink:
			if(a<=0){
				ind=(1:nrowx)>(-a)
				# 'x' shrinked to nothing:
				if(!any(ind)){
					return(NULL)
					}
				else{
					x=x[ind,,drop=FALSE]
					}
				}
			# Extrapolate:
			else{
				x=rbind(outer(a:1,diffx1)+outer(ones(a),x1),x)
				}
			}
		# If 'a' and 'b' specifies the range to which the column given by 'along' is extrapolated/shrinked:
		else{
			toleft=(a-x1[along])/diffx1[along]
			# Shrink:
			if(toleft<1){
				ind=x[,along]>=a
				# 'x' shrinked to nothing:
				if(!any(ind)){
					return(NULL)
					}
				else{
					x=x[ind,,drop=FALSE]
					}
				}
			# Extrapolate:
			else{
				x=rbind(outer(rev(sequence(toleft)),diffx1)+outer(ones(toleft),x1),x)
				}
			}
		nrowx=nrow(x)
		}
	
	# Adding the rest to each end of the matrix:
	if(!grid && nrowx>1 && along!=0){
		restleft=x[1,along]-olda
		if(restleft==0 || identical(numeric(0),restleft) || diffxn[along]!=0){
			restleft=NULL
			}
		retsright=oldb-x[nrowx,along]
		if(retsright==0 || identical(numeric(0),retsright) || diffx1[along]!=0){
			retsright=NULL
			}
		x1=x[1,]
		xn=x[nrowx,]
		diffx1=x1-x[2,]
		diffxn=xn-x[nrowx-1,]
		x=rbind(x1-diffx1/diffx1[along]*restleft,x,xn+diffxn/diffxn[along]*retsright)
		}
	
	
	##### Output #####
	# If byrow==TRUE the columns are mapped back to rows:
	if(byrow){
		x=t(x)
		}
	rownames(x)=NULL
	if(drop.out){
		drop(x)
		}
	else{
		x
		}
	##################################################
	##################################################
	}
