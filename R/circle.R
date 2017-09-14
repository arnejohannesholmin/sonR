#*********************************************
#*********************************************
#' Returns a matrix of two columns representing x- and y-values of positions on a circle of radius 'r' and origin 'origin', defined by the angles given by 'ang' in a polar coordinate system centered at 'origin', or by 'x' or 'y' in a cartesian coordinate system.
#'
#' @param origin  is the origin of the circle.
#' @param r  is the radius of the circle, possibly a vector of the same length as 'ang' for varying radius.
#' @param ang  is a vector of angles defining the point on the circle, as seen in a polar coordinate system centered at 'origin'.
#' @param x  is a vector of x-values at which the circle is caclulated.
#' @param y  is a vector of y-values at which the circle is caclulated.
#' @param n  is (if length(and)==2) the length of 'ang' as a sequence between its first and the second element.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname circle
#'
circle<-function(origin=c(0,0),r=1,ang=seq(0,2*pi,length.out=100),x=NULL,y=NULL,n=100){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-05-16 - Finished.
	# Update:  2009-05-08 - Cleaned up, and translated to english.
	# Last:  2010-08-27 - Replaced data frame output by list output.
	########### DESCRIPTION: ###########
	# Returns a matrix of two columns representing x- and y-values of positions on a circle of radius 'r' and origin 'origin', defined by the angles given by 'ang' in a polar coordinate system centered at 'origin', or by 'x' or 'y' in a cartesian coordinate system.
	########## DEPENDENCIES: ###########
	# contains()
	############ VARIABLES: ############
	# ---origin--- is the origin of the circle.
	# ---r--- is the radius of the circle, possibly a vector of the same length as 'ang' for varying radius.
	# ---ang--- is a vector of angles defining the point on the circle, as seen in a polar coordinate system centered at 'origin'.
	# ---x--- is a vector of x-values at which the circle is caclulated.
	# ---y--- is a vector of y-values at which the circle is caclulated.
	# ---n--- is (if length(and)==2) the length of 'ang' as a sequence between its first and the second element.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Return list if input is a list:
	listinputorigin=FALSE
	# Support for list input:
	if(is.list(origin) && length(origin)>2){
		names(origin)=tolower(names(origin))
		r=origin$r
		ang=origin$ang
		x=origin$x
		y=origin$y
		n=origin$n
		origin=origin$origin
		listinputorigin=TRUE
		}
	
	# 'origin' needs to have length>1:
	origin=rep(origin,length.out=2)
	# If 'x' is given, it overrides 'ang':
	if(!is.null(x)){
		# 'x' and 'r' need to have equal length:
		lx=length(x)
		lr=length(r)
		if(lx!=lr && lr!=1){
			l=min(lx,lr)
			x=x[1:l]
			r=r[1:l]
			}
		# Adding 'origin':
		x=x-origin[1]
		if(max(abs(x))>r){
			x=x[contains(x,c(-r,r),"o")]
			warning("x-values outside of the legal range discarded.")
			}
		y=sqrt(r^2-x^2)
		return(cbind(x+origin[1],y+origin[2]))
		}
	# Else if 'y' is given, it overrides 'ang':
	else if(!is.null(y)){
		# 'y' and 'r' need to have equal length:
		ly=length(y)
		lr=length(r)
		if(ly!=lr && lr!=1){
			l=min(ly,lr)
			y=y[1:l]
			r=r[1:l]
			}
		# Adding 'origin':
		y=y-origin[2]
		if(max(abs(y))>r){
			y=y[contains(y,c(-r,r),"o")]
			warning("y-values outside of the legal range discarded.")
			}
		x=sqrt(r^2-y^2)
		return(cbind(x+origin[1],y+origin[2]))
		}
	# If 'n' is given, the first two elements of ang are used to define a equispaced sequence of length 'n':
	if(length(ang)==2 || n!=100){
		ang=seq(ang[1],ang[length(ang)],length.out=n)
		}
	# 'ang' and 'r' need to have equal length:
	lang=length(ang)
	lr=length(r)
	if(lang!=lr && lr!=1){
		l=min(lang,lr)
		ang=ang[1:l]
		r=r[1:l]
		}
	
	
	##### Execution and output #####
	xy=cbind(origin[1]+r*cos(ang),origin[2]+r*sin(ang))
	colnames(xy)=c("x","y")
	if(listinputorigin){
		list(x=xy[,1],y=xy[,2])
		}
	else{
		xy
		}
	##################################################
	##################################################
	}
