#*********************************************
#*********************************************
#' Returns a function representing the product of three dimensional rotation matrices, specified by the string 'ax'.
#'
#' @param ax  is a string containing a sequence of the characters "x", "y" and "z".
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD repm
#'
#' @export
#' @rdname rotfun
#'
rotfun<-function(ax){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-03-10 - Finished.
	# Update: 2009-03-08 - Generalizing to ignoring other characters in 'ax' than "x", "y" and "z".
	# Last: 2010-02-08 - Changed from using px3x3 matrices to uding 3x3xp matrices and replaced for loops by matrix calculations. Saves time but demands memory!!!
	########### DESCRIPTION: ###########
	# Returns a function representing the product of three dimensional rotation matrices, specified by the string 'ax'.
	########## DEPENDENCIES: ###########
	# repm(), rotx(), roty(), rotz()
	############ VARIABLES: ############
	# ---ax--- is a string containing a sequence of the characters "x", "y" and "z".
	
	
	##################################################
	##################################################
	##### Preparation #####
	# 'charv' is a vector of single characters of the input string 'ax':
	charv=tolower(unlist(strsplit(ax,"")))
	# Removing other characters than "x", "y" and "z".
	xyz=(charv=="x")+(charv=="y")+(charv=="z")
	charv=charv[as.logical(xyz)]
	lcharv=length(charv)
	
	##### Execution #####
	# 'fun' is the function to be returned:
	fun<-function(ang,radians=TRUE,block.out=FALSE,drop.out=TRUE,byrow=TRUE){
		
		########### DESCRIPTION: ###########
		# Returns an array of matrices representing products of three dimensional rotation matrices, as specified by the string 'ax'. (See 'block.out' in the VARIABLES section below.)
		########## DEPENDENCIES: ###########
		# repm(), rotx(), roty(), rotz()
		############ VARIABLES: ############
		# - ('n' is the number of points to be rotated ('lx' in the function rotate()), 'p' is the number of rotations ('lA' in the function rotate() and 'lang' in the function rotfun()), 'q' is the number of axial rotations for each rotation ('lcharv' in the function rotfun()))
		# - ---ax' is a string containing 'q' of the letters 'x', 'y' and 'z', representing the order of the rotations. As an example by="xmy" defines a q=2 xy-rotation (all other letters are discarded).
		# - ---ang' is a vector or matrix of rotation angles. If 'ang' is a vector of length 'lang' and q==1, it is treated as a [lang,1] matrix. If 'ang' is a vector of length 'lang' and q>1, it is transformed in length by the rep(,length.out=q) function to a [1,q] matrix. If 'ang' is a matrix it is, if necessary, transformed in width by the repm() function to a [p,q] matrix.
		# - ---block.out' is TRUE if the output rotation matrices are to be given as a [3,3*p] vertically block divided matrix where the first rotation matrix is A[1:3,1:3], the second is A[1:3,4:6] and so on, and FALSE if the output rotation matrices are to be given as a three dimensional [3,3,p] array where the first rotation matrix is A[,,1], the second is A[,,1] and so on.
		# - ---drop.out' is TRUE if the output rotation matrices 'A' should be stripped of empty dimensions by the drop() function.
		# - ---byrow' is TRUE if the matrices are to be organized by rows (vertically).
		
	
		##### Preparation #####
		# Transformation from degrees to radians:
		if(!radians){
			ang=ang/180*pi
			}
		# If 'ang' is a vector, it needs some treatment using dim() (see description of 'ang' in the VARIABLES section):
		if(length(dim(ang))==0){
			if(lcharv==1){
				dim(ang)=c(length(ang),1)
				}
			else{
				ang=rep(ang,length.out=lcharv)
				dim(ang)=c(1,lcharv)
				}
			}
		# If the number of coloums of 'ang' don't agree with the number of rotations given in 'charv', repm() is used on 'ang'
		dimang=dim(ang)
		if(dimang[2]!=lcharv){
			ang=repm(ang,length.out=lcharv)
			}
		dimang=dim(ang)
		# 'lang' is the same as 'p' in the VARIABLES section:
		lang=dimang[1]
		
		
		##### Execution and output #####
		# Obtaining the output matrix as a product of single rotation matrices:
		# Creating an array of 'lang' square diagonal matrices:
		A=repm(diag(1,3),lang,byrow=TRUE)
		dimA=dim(A)
		if(lcharv==0){
			return(A)
			}
		# For each letter in charv we multiply with the result from rotx(), roty() or rotz():
		for(j in 1:lcharv){
			if(identical(charv[j],"x")){
				A0=repm(rotx(ang[,j],block.out=TRUE),times=3,byrow=TRUE)
				A=rep(A,each=3)
				dim(A)=c(3,9*lang)
				A=blockshift(A,c(3,3))
				A=A0 * A
				A=rowSums(A)
				}
			if(identical(charv[j],"y")){
				A0=repm(roty(ang[,j],block.out=TRUE),times=3,byrow=TRUE)
				A=rep(A,each=3)
				dim(A)=c(3,9*lang)
				A=blockshift(A,c(3,3))
				A=A0 * A
				A=rowSums(A)
				}
			if(identical(charv[j],"z")){
				A0=repm(rotz(ang[,j],block.out=TRUE),times=3,byrow=TRUE)
				A=rep(A,each=3)
				dim(A)=c(3,9*lang)
				A=blockshift(A,c(3,3))
				A=A0 * A
				A=rowSums(A)
				}
			}
		dim(A)=dimA
			
		if(block.out){
			if(!byrow){
				A=blockshift(A,c(3,3))
				}
			}
			
		if(!block.out){
			if(drop.out && lang!=1){
				dim(A)=c(3,lang,3)
				A=aperm(A,perm=c(1,3,2))
				}
			}
		A
		}
	
	
	##### Output #####
	fun
	##################################################
	##################################################
	}
