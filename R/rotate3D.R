#*********************************************
#*********************************************
#' Rotates the set of 'n' three dimensional points 'x' by the 'p' rotation matrices given by 'A', or by the rotation angles 'ang' around the axes given by 'by'. More than one rotations may be given, either as a [3,3,p] array holding 'p' rotation matrices, a [3,3*p] vertically block divided matrix also holding 'p' rotation matrices, or as 'p' pairs of significant characters "x", "y" and "z" and rotation angles given by 'ang'. When p>1, the output has dimensions [n,3,p], where 'n' is the number of points to be rotated. When paired==TRUE, the output has dimension [n,3].
#'
#' @param x  is the input given as a matrix of columns representing the x-, y- and z-coordinates (size: [n,3]), or a matrix of rows representing the x-, y- and z-coordinates (size: [3,n]) in which case 'paired' needs to be TRUE, or a list of elements named "x", "y" and "z".
#' @param by  is a string containing 'q' of the letters 'x', 'y' and 'z', representing the order of the rotations. As an example by="x and y" defines a q=2 xy-rotation (all other letters are discarded).
#' @param ang  is a vector or matrix of rotation angles. If 'ang' is a vector of length 'lang' and q==1, it is treated as a [lang,1] matrix. If 'ang' is a vector of length 'lang' and q>1, it is transformed in length by the rep(,length.out=q) function to a [1,q] matrix. If 'ang' is a matrix it is, if necessary, transformed in width by the repm() function to a [p,q] matrix.
#' @param paired  is set to TRUE if each point given by 'x' should be rotated by the corresponding rotation matrix defined by 'A' or by 'by' and 'ang'. If n!=p, the shorter is expanded to the longer. 
#' @param radians  is TRUE if the input angles 'ang' are in radians, and FALSE if degrees are used.
#' @param sph.in  is TRUE if the input points 'x' are in spherical coordinates.
#' @param sph.out  is TRUE if the output points should be in spherical coordinates.
#' @param drop.out  is TRUE if the output should be stripped of empty dimensions by the drop() function.
#' @param list.out  is FALSE to force array output, TRUE to force list output and NULL to let the nature of the input decide (list giving list output).
#' @param reverse  is TRUE to reverse the rotation by reversing the order of 'by' and 'ang', and also multiplying 'ang' by -1.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD car2sph repm sph2car zeros
#'
#' @export
#' @rdname rotate3D
#'
rotate3D<-function(x, by, ang, paired=FALSE, radians=TRUE, sph.in=FALSE, sph.out=FALSE, drop.out=TRUE, list.out=NULL, reverse=FALSE, na0=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2014-10-12 - New function using c++, replacing rotate().
	########### DESCRIPTION: ###########
	# Rotates the set of 'n' three dimensional points 'x' by the 'p' rotation matrices given by 'A', or by the rotation angles 'ang' around the axes given by 'by'. More than one rotations may be given, either as a [3,3,p] array holding 'p' rotation matrices, a [3,3*p] vertically block divided matrix also holding 'p' rotation matrices, or as 'p' pairs of significant characters "x", "y" and "z" and rotation angles given by 'ang'. When p>1, the output has dimensions [n,3,p], where 'n' is the number of points to be rotated. When paired==TRUE, the output has dimension [n,3].
	########## DEPENDENCIES: ###########
	# dim_all(), sph2car(), rotfun(), repm(), zeros(), car2sph()
	############ VARIABLES: ############
	# - ('n' is the number of points to be rotated ('lx' in the function rotate()), 'p' is the number of rotations ('lA' in the function rotate() and 'lang' in the function rotfun()), 'q' is the number of axial rotations for each rotation ('lcharv' in the function rotfun()))
	#
	# ---x--- is the input given as a matrix of columns representing the x-, y- and z-coordinates (size: [n,3]), or a matrix of rows representing the x-, y- and z-coordinates (size: [3,n]) in which case 'paired' needs to be TRUE, or a list of elements named "x", "y" and "z".
	# ---by--- is a string containing 'q' of the letters 'x', 'y' and 'z', representing the order of the rotations. As an example by="x and y" defines a q=2 xy-rotation (all other letters are discarded).
	# ---ang--- is a vector or matrix of rotation angles. If 'ang' is a vector of length 'lang' and q==1, it is treated as a [lang,1] matrix. If 'ang' is a vector of length 'lang' and q>1, it is transformed in length by the rep(,length.out=q) function to a [1,q] matrix. If 'ang' is a matrix it is, if necessary, transformed in width by the repm() function to a [p,q] matrix.
	# ---paired--- is set to TRUE if each point given by 'x' should be rotated by the corresponding rotation matrix defined by 'A' or by 'by' and 'ang'. If n!=p, the shorter is expanded to the longer. 
	# ---radians--- is TRUE if the input angles 'ang' are in radians, and FALSE if degrees are used.
	# ---sph.in--- is TRUE if the input points 'x' are in spherical coordinates.
	# ---sph.out--- is TRUE if the output points should be in spherical coordinates.
	# ---drop.out--- is TRUE if the output should be stripped of empty dimensions by the drop() function.
	# ---list.out--- is FALSE to force array output, TRUE to force list output and NULL to let the nature of the input decide (list giving list output).
	# ---reverse--- is TRUE to reverse the rotation by reversing the order of 'by' and 'ang', and also multiplying 'ang' by -1.
	
	
	##################################################
	##################################################
	##### Preparation #####
	###The max amount of memory occupied by the function:  ###
	# Return NULL is the input is NULL:
	if(length(x)==0){
		return()
		}
	
	# If the input points are to be transformed from spherical to cartesian coordinates:
	if(sph.in){
		x <- sph2car(x)
		}
		
	# Support for list input for 'x':
	listinput <- FALSE
	if(is.list(x)){
		listinput <- TRUE
		names(x) <- tolower(names(x))
		if(length(x$x)>0 && length(x$y)>0 && length(x$z)>0){
			x <- cbind(x$x,x$y,x$z)
		}
		else{
			x <- cbind(x[[1]],x[[2]],x[[3]])
		}
	}
	# Array input for 'x':
	else if(length(dim(x))==2){
		if(ncol(x)>3){
			warning("Only the first 3 columns of x used.")
			x <- x[,1:3]
		}
		else if(ncol(x)<3){
			warning("The matrix x padded by NAs to a 3 columns matrix.")
			x <- cbind(x, array(NA, dim=c(1, 3-ncol(x))))
		}
	}
	else{
		if(length(x)<3){
			stop("Three coordinates must be given")
		}
		x <- matrix(x[1:3],nrow=1,ncol=3)
	}
	
	
	##### Execution #####
	if(!is.double(x)){
		dimx <- dim(x)
		x <- as.double(x)
		dim(x) <- dimx
		}
	Npoints <- nrow(x)
	
	charv <- tolower(unlist(strsplit(by,"")))
	# Removing other characters than "x", "y" and "z".
	xyz <- (charv=="x")+(charv=="y")+(charv=="z")
	charv <- charv[as.logical(xyz)]
	by <- as.integer(match(charv,c("x","y","z"))-1)
	Lby <- length(by)
	if(reverse){
		by <- rev(by)
		}
	
	# If 'ang' is a vector, it needs some treatment using dim() (see description of 'ang' in the VARIABLES section):
	dimang <- dim(ang)
	if(length(dimang)==0){
		if(Lby==1){
			dim(ang) <- c(length(ang),1)
			}
		else{
			ang <- rep(ang,length.out=Lby)
			dim(ang) <- c(1,Lby)
			}
		}
	dimang <- dim(ang)
	# If the number of coloums of 'ang' don't agree with the number of rotations given in 'charv', repm() is used on 'ang'
	if(dimang[2]!=Lby){
		ang <- repm(ang,length.out=Lby)
		}
	if(!is.double(ang)){
		dimang <- dim(ang)
		ang <- as.double(ang)
		dim(ang) <- dimang
		}
	Nrot <- NROW(ang)
	# Test for missing values
	if(na0){
		ang[is.na(ang)] <- 0
		}
	
	if(reverse){
		if(length(dim(ang))==2){
			ang <- -ang[,ncol(ang):1]
			}
		else{
			ang <- -rev(ang)
			}
		}
	
	A <- zeros(3,3,Nrot)
	A[1,1,] <- 1
	A[2,2,] <- 1
	A[3,3,] <- 1
	
	
	# Rotate using c++:
	# Allow for NAs:
	isNAx <- is.na(x)
	x[isNAx] <- 0
	U <- .C("rotate3d", x, Npoints, by, Lby, ang, A, Nrot, as.integer(paired), rep(x,if(paired) 1 else Nrot), PACKAGE="sonR")
	U[[9]][isNAx] <- NA
	dim(U[[6]]) <- dim(A)
	if(paired){
		dim(U[[9]]) <- c(Npoints,3)
		}
	else{
		dim(U[[9]]) <- c(Npoints,3,Nrot)
		}
	
	# If 'paired' is FALSE an array of three dimensions is returned:
	if(sph.out){
		if(paired){
			U[[9]] <- car2sph(U[[9]])
			}
		else{
			tempout <- car2sph(c(U[[9]][,1,]),c(U[[9]][,2,]),c(U[[9]][,3,]))
			U[[9]][,1,] <- tempout[,1]
			U[[9]][,2,] <- tempout[,2]
			U[[9]][,3,] <- tempout[,3]
			}
		}
	
	
	##### Output #####
	# If input 'x' was a list, a list is returned:
	if(isTRUE(list.out) || (listinput && !identical(list.out,FALSE))){
		if(length(dim(U[[9]]))==3){
			U[[9]] <- list(x=U[[9]][,1,,drop=drop.out],y=U[[9]][,2,,drop=drop.out],z=U[[9]][,3,,drop=drop.out])
			}
		else{
			U[[9]] <- list(x=U[[9]][,1,drop=drop.out],y=U[[9]][,2,drop=drop.out],z=U[[9]][,3,drop=drop.out])
			}
		}
	if(drop.out){
		drop(U[[9]])
		}
	else{
		U[[9]]
		}
	##################################################
	##################################################
	}
