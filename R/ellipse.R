#*********************************************
#*********************************************
#' Returns a matrix of two columns representing x- and y-values of positions on an ellipse of horizontal axis 'a', vertical axis 'b' and origin 'origin', defined by the angles given by 'ang' in a polar coordinate system centered at 'origin', or by 'x' or 'y' in a cartesian coordinate system.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname ellipse
#'
ellipse <- function(origin=c(0,0), a=2, b=1, ang=seq(0, 2*pi, length.out=n), x=NULL, y=NULL, n=100, aspectcomp=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-05-16 - Finished.
	# Update:  2009-05-08 - Cleaned up, and translated to english.
	# Update:  2009-08-02 - Added support for varying 'a' and 'b' values.
	# Last:  2010-08-27 - Replaced data frame output by list output.
	########### DESCRIPTION: ###########
	# Returns a matrix of two columns representing x- and y-values of positions on an ellipse of horizontal axis 'a', vertical axis 'b' and origin 'origin', defined by the angles given by 'ang' in a polar coordinate system centered at 'origin', or by 'x' or 'y' in a cartesian coordinate system.
	########## DEPENDENCIES: ###########
	# contains()
	############ VARIABLES: ############
	# - 'origin' is the origin of the ellipse.
	# - 'a' is the horizontal axis of the ellipse.
	# - 'b' is the vertical axis of the ellipse.
	# - 'ang' is a vector of angles defining the points on the ellipse, as seen in a polar coordinate system centered at 'origin'.
	# - 'x' is a vector of x-values at which the ellipse is caclulated.
	# - 'y' is a vector of y-values at which the ellipse is caclulated.
	# - 'n' is (if different from NULL) the length of 'ang' as a sequence between the first and the second element of 'ang'.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Return list if input is a list:
	listinputorigin = FALSE
	# Support for list input:
	if(is.list(origin) && length(origin)>2){
		names(origin) = tolower(names(origin))
		a = origin$a
		b = origin$b
		ang = origin$ang
		x = origin$x
		y = origin$y
		n = origin$n
		origin = origin$origin
		listinputorigin = TRUE
		}
	
	# 'origin' needs to have length>1:
	origin = rep(origin, length.out=2)
	# If 'x' is given, it overrides 'ang':
	if(!is.null(x)){
		# 'x', 'a' and 'b' need to have equal length:
		lx = length(x)
		la = length(a)
		lb = length(b)
		if((lx!=la || lx!=lb) && la!=1 && lb!=1){
			l = min(lx, la, lb)
			x = x[1:l]
			a = a[1:l]
			b = b[1:l]
			}
		# Adding 'origin':
		x = x-origin[1]
		if(max(abs(x))>a){
			x = x[contains(x, c(-a,a), "o")]
			warning("x-values outside of the legal range are discarded.")
			}
		y = sqrt(b^2*(1-(x/a)^2))
		return(cbind(x+origin[1], y+origin[2]))
		}
	# Else if 'y' is given, it overrides 'ang':
	else if(!is.null(y)){
		# 'y', 'a' and 'b' need to have equal length:
		ly = length(y)
		la = length(a)
		lb = length(b)
		if((ly!=la || ly!=lb) && la!=1 && lb!=1){
			l = min(ly, la, lb)
			y = y[1:l]
			a = a[1:l]
			b = b[1:l]
			}
		# Adding 'origin':
		y = y-origin[2]
		if(max(abs(y))>b){
			y = y[contains(y, c(-b,b), "o")]
			warning("y-values outside of the legal range are discarded.")
			}
		x = sqrt(a^2*(1-(y/b)^2))
		return(cbind(x+origin[1], y+origin[2]))
		}
	# If 'n' is given, the first two elements of ang are used to define a equispaced sequence of length 'n':
	#if(length(ang)==2 || n!=100){
	#	ang = seq(ang[1], ang[length(ang)], length.out=n)
	#	}
	if(length(ang)==2){
		ang = seq(ang[1], ang[2], length.out=n)
	}
	if(aspectcomp){
		x <- seq(0, 2*pi, l=n)
		if(a>b){
			diffang <- abs(sin(x))^(0.05 * a/b)
		}
		else{
			diffang <- abs(cos(x))^(0.05 * b/a)
		}
		ang <- cumsum(diffang) / sum(diffang) * 2*pi
	}	
		
		
	# 'ang', 'a' and 'b' need to have equal length:
	lang = length(ang)
	la = length(a)
	lb = length(b)
	if((lang!=la || lang!=lb) && la!=1 && lb!=1){
		l = min(lang, la, lb)
		ang = ang[1:l]
		a = a[1:l]
		b = b[1:l]
		}
	
	
	##### Execution and output #####
	eps2 = (a^2-b^2)/a^2
	rho = sqrt(b^2/ (1 - eps2 * cos(ang)^2) )
	xy = cbind(origin[1] + rho * cos(ang), origin[2] + rho * sin(ang))
	colnames(xy) = c("x","y")
	if(listinputorigin){
		list(x=xy[,1], y=xy[,2])
		}
	else{
		xy
		}
	##################################################
	##################################################
	}
