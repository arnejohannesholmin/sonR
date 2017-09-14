#*********************************************
#*********************************************
#' Returns three dimensional rotation matrices for rotation around the z-axis of a three dimensional coordinate system.
#'
#' @param ang  is a vedctor of the rotation angles around the z-axis.
#' @param radians  is TRUE if 'ang' is given in radians, and FALSE if 'ang' is given in degrees.
#' @param block.out  is TRUE if the rotation matrices should be returned as a 3x3*l matrix, where 'l' is the number of rotation angles. If FALSE, an 3x3xl array is returned.
#' @param drop.out  is relevant only if block.out==FALSE and l==1, and if set to TRUE a 3x3 matrix is returned instead of a 3x3x1 array.
#' @param byrow  is TRUE if the matrices are to be organized by rows (vertically).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname rotz
#'
rotz<-function(ang,radians=TRUE,block.out=FALSE,drop.out=TRUE,byrow=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-03-10 - Finished.
	# Update: 2009-03-08 - Simplified version.
	# Update: 2010-02-08 - Changed from using px3x3 matrices to uding 3x3xp matrices.
	# Last: 2010-06-10 - Removed dependency on blockshift().
	########### DESCRIPTION: ###########
	# Returns three dimensional rotation matrices for rotation around the z-axis of a three dimensional coordinate system.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---ang--- is a vedctor of the rotation angles around the z-axis.
	# ---radians--- is TRUE if 'ang' is given in radians, and FALSE if 'ang' is given in degrees.
	# ---block.out--- is TRUE if the rotation matrices should be returned as a 3x3*l matrix, where 'l' is the number of rotation angles. If FALSE, an 3x3xl array is returned.
	# ---drop.out--- is relevant only if block.out==FALSE and l==1, and if set to TRUE a 3x3 matrix is returned instead of a 3x3x1 array.
	# ---byrow--- is TRUE if the matrices are to be organized by rows (vertically).
	
	
	##################################################
	##################################################
	if(!radians){
		ang=ang/180*pi
		}
	l=length(ang)
	# Mehtods for byrow=TRUE and byrow=FALSE differ only be the minus sign on the sine terms and the transpose
	if(block.out){
		if(byrow){
			a=rbind(cos(ang),sin(ang),0,-sin(ang),cos(ang),0,0,0,1)
			dim(a)=c(3,3*l)
			a=t(a)
			}
		else{
			a=rbind(cos(ang),-sin(ang),0,sin(ang),cos(ang),0,0,0,1)
			dim(a)=c(3,3*l)
			}
		}
	else{
		a=rbind(cos(ang),-sin(ang),0,sin(ang),cos(ang),0,0,0,1)
		dim(a)=c(3,3,l)
		}
	# Output:
	if(drop.out){
		drop(a)
		}
	else{
		a
		}
	##################################################
	##################################################
	}
