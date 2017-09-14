#*********************************************
#*********************************************
#' Returns a list of x-values and y-values of points defined by the fractions 'w' along the partial circular path defined by the list or matrix of origin points 'origin', the vector of radii 'r' and the list or matrix of angles 'ang'. Cumulative lengths of the path are also returned.
#'
#' @param origin  is the origin of the path, given as a list of two elements or a matrix of two columns, representing the x and y values of the origins.
#' @param r  is a vector of length 'n' holding the input radius-variable.
#' @param ang  is the input angle-variable as defined on the unit circle, given as a list of two elements or a matrix of two columns, representing the start and the end angles of the circle segments.
#' @param w  are the positions along the path (x,y), or along the x-axis if alongx==TRUE. If 'w' is a single integer, it is set to seq(0,1,length.out=w). If 'speed' is given 'w' is interpreted as time values.
#' @param dw  are the piecewise lengths along the path (x,y), or along the x-axis if alongx==TRUE, used only if w!=NULL. If 'dw' ihas length>1 it overrides 'w' by w=cumsum(dw).
#' @param speed  is a vector of sound speeds of each circle segment. If given, 'w' is interpreted as time values starting from dt=0 at the start of the path.
#' @param normalized  is TRUE if the positions given by 'w' (and 'speed') are or should be normalized to the set [0,1], where 1 represents the end of the path.
#' @param list.out  is TRUE if the output should be a data frame with the names "x", "y" and "l". Else a matrix with the same names is returned.
#' @param plot  is TRUE if a simple plot is to be drawn.
#' @param ...  are parameters to be passed on to plot() if plot==TRUE.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD setrange zeros
#' @importFrom utils tail
#' @importFrom graphics hist lines points
#'
#' @export
#' @rdname track.circle
#'
track.circle<-function(origin,r,ang,w=NULL,dw=NULL,speed=NULL,normalized=FALSE,list.out=TRUE,plot=FALSE,...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-06-12 - First clean version.
	# Update: 2010-02-19 - Added suport for time differences 'dw' if speed is given.
	# Last:  2010-08-26 - Removed the data.frame output and expanded to return 'pos' as is done in track.line().
	########### DESCRIPTION: ###########
	# Returns a list of x-values and y-values of points defined by the fractions 'w' along the partial circular path defined by the list or matrix of origin points 'origin', the vector of radii 'r' and the list or matrix of angles 'ang'. Cumulative lengths of the path are also returned.
	########## DEPENDENCIES: ###########
	# setrange(), zeros(), circle()
	############ VARIABLES: ############
	# ('n' is the number of circle segments)
	# ---origin--- is the origin of the path, given as a list of two elements or a matrix of two columns, representing the x and y values of the origins.
	# ---r--- is a vector of length 'n' holding the input radius-variable.
	# ---ang--- is the input angle-variable as defined on the unit circle, given as a list of two elements or a matrix of two columns, representing the start and the end angles of the circle segments.
	# ---w--- are the positions along the path (x,y), or along the x-axis if alongx==TRUE. If 'w' is a single integer, it is set to seq(0,1,length.out=w). If 'speed' is given 'w' is interpreted as time values.
	# ---dw--- are the piecewise lengths along the path (x,y), or along the x-axis if alongx==TRUE, used only if w!=NULL. If 'dw' ihas length>1 it overrides 'w' by w=cumsum(dw).
	# ---speed--- is a vector of sound speeds of each circle segment. If given, 'w' is interpreted as time values starting from dt=0 at the start of the path.
	# ---normalized--- is TRUE if the positions given by 'w' (and 'speed') are or should be normalized to the set [0,1], where 1 represents the end of the path.
	# ---list.out--- is TRUE if the output should be a data frame with the names "x", "y" and "l". Else a matrix with the same names is returned.
	# ---plot--- is TRUE if a simple plot is to be drawn.
	# ---...--- are parameters to be passed on to plot() if plot==TRUE.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Support for list input and matrix input for 'origin':
	if(is.list(origin) && length(origin)>1){
		names(origin)=tolower(names(origin))
		if(!is.null(origin$x) && !is.null(origin$y)){
			ox=origin$x
			oy=origin$y
			}
		else{
			ox=origin[[1]]
			oy=origin[[2]]
			}
		}
	else if(length(dim(origin))==2){
		ox=origin[,1]
		oy=origin[,2]
		}
	
	# Support for list input and matrix input for 'ang':
	if(is.list(ang) && length(ang)>1){
		names(ang)=tolower(names(origin))
		if(!is.null(ang$a) && !is.null(ang$b)){
			ang1=ang$a
			ang2=ang$b
			}
		else{
			ang1=ang[[1]]
			ang2=ang[[2]]
			}
		}
	else if(length(dim(ang))==2){
		ang1=ang[,1]
		ang2=ang[,2]
		}
		
	# 'r' must be non-negative:
	if(min(r)<0){
		warning("'r' must be non-negative (abs(r) used)")
		r=abs(r)
		}	
	
	# Lengths of the inputs. Must be equal:
	lx=length(ox)
	lr=length(r)
	lang=length(ang1)
	if(!all(c(lx,lr,lang)==lx)){
		stop("Input lengths must agree. See info(\"track.circle\")")
		}
	
	
	##### Execution #####
	# 'angdiff' and 'leftright' are the span of the circle segments and the direction on the circle segment defined as clockwise = -1 and counter clockwise = 1.
	angdiff=ang2-ang1
	leftright=sign(angdiff)
	# The length variable to be used in the function:
	lengths=abs(r*angdiff)
	cumlengths=c(0,cumsum(lengths))
	totlengths=cumlengths[lx+1]
	
	# If the length of 'dw' exceeds 1, it overrides 'w':
	if(length(dw)>1){
		w=cumsum(dw)
		}
	# If w==NULL, regularly spaced points along the path are returned:
	else if(is.null(w)){
		w=seq(0,1,length.out=lx)
		normalized=TRUE
		}
	# If 'w' has length 1 regularly spaced points are returned:
	else if(length(w)==1){
		if(!is.null(dw)){
			w=seq(0,by=dw,length.out=w)
			}
		else{
			w=seq(0,1,length.out=w)
			normalized=TRUE
			}
		}
	# 'w' needs to be sorted:
	else{
		w=sort(w)
		}
	
	# If 'speed' is given, the calculation is done on times not lengths:
	if(!is.null(speed)){
		t=w
		# The length variable to be used in the function:
		times=abs(r*angdiff)/speed
		cumtimes=c(0,cumsum(times))
		tottimes=tail(cumtimes,1)
		# Scaling 'w' to match 'cumlengths':
		if(normalized){
			t=setrange(t)*tottimes
			}
		# 'w' must be non-negative:
		if(min(t)<0 || max(t)>tottimes){
			warning("Elements of 'w' (time) exceed the path defined by the given circles, and are ignored")
			t=t[t>=0 & t<=tottimes]
			}
		
		# 'pos' provides in which line segments the points defined by 'w' are positioned. If 'totheend' is TRUE, special attention is needed in the calculation to avoid an error:
		pos=hist(t,breaks=cumtimes,plot=FALSE)$counts
		pos=rep(1:lx,pos)
		# Length of the circle segments to the left of the points:
		leftang=(t-cumtimes[pos])*speed[pos]/(r[pos])
		
		# The output consists of x-positions 'x', y-positions 'y', lengths 'l' and time/length 'w' along the circular path:
		w=cumlengths[pos]+r[pos]*leftang
		xywout=cbind(ox[pos]+r[pos]*cos(ang1[pos]+leftang*leftright[pos]),oy[pos]+r[pos]*sin(ang1[pos]+leftang*leftright[pos]),w,t,pos)
		colnames(xywout)=c("x","y","l","t","pos")
		}
	else{
		# Scaling 'w' to match 'cumlengths':
		if(normalized){
			w=setrange(w)*totlengths
			}
		# 'w' must be non-negative:
		if(min(w)<0 || max(w)>totlengths){
			warning("Elements of 'w' exceed the path defined by the given circles, and are ignored")
			w=w[w>=0 & w<=totlengths]
			}
		# 'pos' provides in which line segments the points defined by 'w' are positioned. If 'totheend' is TRUE, special attention is needed in the calculation to avoid an error:
		pos=hist(w,breaks=cumlengths,plot=FALSE)$counts
		pos=rep(1:lx,pos)
		# Length of the circle segments to the left of the points:
		leftang=(w-cumlengths[pos])/(r[pos])
		
		# The output consists of x-positions 'x', y-positions 'y', lengths 'l' and time/length 'w' along the circulat path:
		xywout=cbind(ox[pos]+r[pos]*cos(ang1[pos]+leftang*leftright[pos]),oy[pos]+r[pos]*sin(ang1[pos]+leftang*leftright[pos]),w,pos)
		colnames(xywout)=c("x","y","l","pos")
		}
		
		
	##### Output #####
	# Plot of the circle segments and the points:
	if(plot){
		dens=100
		circles=zeros(lx,dens,2)
		for(i in 1:lx){
			circles[i,,]=circle(c(ox[i],oy[i]),r[i],seq(ang1[i],ang2[i],length.out=dens))
			}
		plot(NULL,xlim=range(xywout[,1],circles[,,1]),ylim=range(xywout[,2],circles[,,2]),...)
		
		for(i in 1:lx){
			lines(circles[i,,])
			}
		points(xywout,col=2,pch="*",cex=1)
		}
		
	# The output:
	if(list.out){
		if(ncol(xywout)==4){
			list(x=xywout[,1],y=xywout[,2],l=xywout[,3],pos=xywout[,4])
			}
		else if(ncol(xywout)==5){
			list(x=xywout[,1],y=xywout[,2],l=xywout[,3],t=xywout[,4],pos=xywout[,5])
			}
		}
	else{
		xywout
		}
	##################################################
	##################################################
	}
