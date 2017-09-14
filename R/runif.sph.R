#*********************************************
#*********************************************
#' Generates uniformly distributed variables in the spherical volume segments defined by 'r', 'theta' and 'phi'.
#'
#' @param n  is the number of positions to draw, given either as a single integer or as a vector where each element corresponds to a spherical segment given by the rows of 'r', 'theta' and 'phi'.
#' @param r  is range of radial coordinates, given as a vector of interval values (in which case length(r)-1 spherical segemnts are treated) or as a two column matrix where the first column represents the lower limit of the sperical segments, and the second column represents the upper limit.
#' @param theta  is the range of azimuth angles in radians, given in the same way as 'r'.
#' @param phi  is the range of elevation angles in radians, given in the same way as 'r'.
#' @param offset  is the angluar offset by which the input specifications 'r', 'theta' and 'phi' are rotated prior to application of runif.sph(). The result of runif.sph() is rotated back by the inverse angles.
#' @param scale.voxels  is a factor by which to enlarge the extent of each spherical element inside which points are drawn. If scale.voxels=2, the radial positions are drawn from -0.5 to 1.5 of the extent of the radial boundaries, so that point extend into the previous and next radial interval. The same applies to all three directions.
#' @param rand.gen  is either a string naming the random generating function to use when drawing positions inside the spherical segments, or a (symetrical) function. If rand.gen[1]=="beta", the symmetrical beta function is used (shape1=shape2=2).
#' @param reverse  is TRUE to reverse the rotation by reversing the order of 'by' and 'ang', and also multiplying 'ang' by -1.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD repm sph2car strff
#' @importFrom stats rbeta rnorm runif
#'
#' @export
#' @rdname runif.sph
#'
runif.sph<-function(n, r=c(0,1), theta=c(0,2*pi), phi=c(0,pi), sph.out=TRUE, offset=list(), scale.voxels=1, rand.gen=c("unif","beta","norm"), reverse=FALSE){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-05-08 - Clean version.
	# Update: 2011-09-17 - Cleaned up documentation and changed to accept all legnths of 'n' (repeated to the desired length) and to repeat 'r', 'theta' and 'phi' to have the same number of rows.
	# Update: 2013-09-12 - Added the parameter 'offset'.
	# Update: 2013-09-12 - Added the parameters 'scale.voxels' and 'rand.gen'.
	# Update: 2014-11-12 - Changed to use rotate3D().
	# Last: 2015-04-30 - Added support for angles outside of the range [0,2*pi].
	########### DESCRIPTION: ###########
	# Generates uniformly distributed variables in the spherical volume segments defined by 'r', 'theta' and 'phi'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---n--- is the number of positions to draw, given either as a single integer or as a vector where each element corresponds to a spherical segment given by the rows of 'r', 'theta' and 'phi'.
	# ---r--- is range of radial coordinates, given as a vector of interval values (in which case length(r)-1 spherical segemnts are treated) or as a two column matrix where the first column represents the lower limit of the sperical segments, and the second column represents the upper limit.
	# ---theta--- is the range of azimuth angles in radians, given in the same way as 'r'.
	# ---phi--- is the range of elevation angles in radians, given in the same way as 'r'.
	# ---offset--- is the angluar offset by which the input specifications 'r', 'theta' and 'phi' are rotated prior to application of runif.sph(). The result of runif.sph() is rotated back by the inverse angles.
	# ---scale.voxels--- is a factor by which to enlarge the extent of each spherical element inside which points are drawn. If scale.voxels=2, the radial positions are drawn from -0.5 to 1.5 of the extent of the radial boundaries, so that point extend into the previous and next radial interval. The same applies to all three directions.
	# ---rand.gen--- is either a string naming the random generating function to use when drawing positions inside the spherical segments, or a (symetrical) function. If rand.gen[1]=="beta", the symmetrical beta function is used (shape1=shape2=2).
	# ---reverse--- is TRUE to reverse the rotation by reversing the order of 'by' and 'ang', and also multiplying 'ang' by -1.
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(length(scale.voxels)<3){
		scale.voxels = rep(scale.voxels, l=3)
		}
	
	if(!is.function(rand.gen)){
		if(strff("unif",rand.gen[[1]])){
			rand.gen=runif
			}
		else if(strff("beta",rand.gen[[1]])){
			rand.gen=function(x) rbeta(x, 2, 2)
			}
		else if(strff("norm",rand.gen[[1]])){
			rand.gen=function(x) rnorm(x, 0.5)
			}
		}
	
	# Issue a warning if any r<0, theta<0, theta>2*pi, phi<0, phi>pi:
	if(length(r) && min(r,na.rm=TRUE)<0){
		warning("Negative values of 'r', possibly leading to negative r-values in the output")
		}
	
	# 'n' must be an integer:
	n=floor(n)
	
	# Accept if only one value is given:
	if(length(r)==1){
		r=c(r,r)
		}
	if(length(theta)==1){
		theta=c(theta,theta)
		}
	if(length(phi)==1){
		phi=c(phi,phi)
		}
	# If given as vectors, transform to matrices of intervals:
	if(length(dim(r))==0){
		r=cbind(r[-length(r)],r[-1])
		}
	if(length(dim(theta))==0){
		theta=cbind(theta[-length(theta)],theta[-1])
		}
	if(length(dim(phi))==0){
		phi=cbind(phi[-length(phi)],phi[-1])
		}
	
	# Accept only angles on the unit sphere:
	#if(any(theta<0 | theta>2*pi)){
	#	warning(paste("The azimuth angle 'theta' must be in [0,2*pi]. Values outside of this range were set to the closest range value. Range of theta:",paste(range(theta),collapse=", ")))
	#	theta[theta<0]=0
	#	theta[theta>2*pi]=2*pi
	#	}
	#if(any(phi<0 | phi>pi)){
	#	warning(paste("The elevation angle 'phi' must be in [0,pi]. Values outside of this range were set to the closest range value. Range of phi:",paste(range(phi),collapse=", ")))
	#	phi[phi<0]=0
	#	phi[phi>pi]=pi
	#	}
		
	# Ensure that 'r', 'theta' and 'phi' have equal dimension:
	if(diff(range(c(nrow(r),nrow(theta),nrow(phi)))) <= 0){
		m=max(c(nrow(r),nrow(theta),nrow(phi)))
		r=repm(r,length.out=m,byrow=TRUE)
		theta=repm(theta,length.out=m,byrow=TRUE)
		phi=repm(phi,length.out=m,byrow=TRUE)
		}
	# If 'n' does not have length equal to the number of rows of 'r', repeat to the desired length:
	if(length(n)!=nrow(r)){
		n=rep(n,length.out=nrow(r))
		}
	
	N=sum(n, na.rm=TRUE)
	# Return NULL if N==0:
	if(N==0){
		return()
		}
	# Repeat the values of 'r', 'theta' and 'phi' for use in the generation of random positions:
	#r=cbind(rep(r[,1],n),rep(r[,2],n))
	#theta=cbind(rep(theta[,1],n),rep(theta[,2],n))
	#phi=cbind(rep(phi[,1],n),rep(phi[,2],n))
	# Draw the random bases for the spherical random positions:
	#ur=setrange(rand.gen(N),0.5+c(-scale.voxels/2,scale.voxels/2),check.ends=FALSE)
	#utheta=setrange(rand.gen(N),0.5+c(-scale.voxels/2,scale.voxels/2),check.ends=FALSE)
	#uphi=setrange(rand.gen(N),0.5+c(-scale.voxels/2,scale.voxels/2),check.ends=FALSE)
	# THE ABOVE CODE WAS A BUG, AS WE DO NOT WANT THE MINIMUM AND MAXIMUM TO BE EXACTLY AT scale.voxels*c(0,1), BUT SIMPLY SCALED AROUND THE MIDPOINT 0.5:
	ur=rand.gen(N)
	ur=scale.voxels[1]*(ur-0.5)+0.5
	utheta=rand.gen(N)
	utheta=scale.voxels[1]*(utheta-0.5)+0.5
	uphi=rand.gen(N)
	uphi=scale.voxels[1]*(uphi-0.5)+0.5
	
	
	##### Execution #####
	# Apply the randomly generated values 'ur' to the radial dimension:
	r1=rep(r[,1],n)
	r1 = r1*r1*r1
	r2=rep(r[,2],n)
	r2 = r2*r2*r2
	r=( r1 + ur*(r2 - r1) )^(1/3)
	rm(r1)
	rm(r2)
	rm(ur)
	
	# Apply the randomly generated values 'utheta' to the azimuth dimension:
	theta=rep(theta[,1],n) + utheta*(rep(theta[,2],n)-rep(theta[,1],n))
	rm(utheta)
	
	# Apply the randomly generated values 'uphi' to the elevation dimension:
	cosphi1=cos(rep(phi[,1],n))
	cosphi2=cos(rep(phi[,2],n))
	phi=acos(cosphi1-uphi*(cosphi1 - cosphi2))
	
	#phi=acos(cos(rep(phi[,1],n)) - uphi*( cos(rep(phi[,1],n)) - cos(rep(phi[,2],n)) ))
	rm(uphi)
	rm(cosphi1)
	rm(cosphi2)
	
	
	##### Output #####
	# If 'offset' is given, rotate back to the proper orientation:
	if(length(offset)>0){
		#rotate3D(cbind(r=r,theta=theta,phi=phi), by=rev(offset$by), ang=-offset$ang, sph.in=TRUE, sph.out=sph.out)
		rotate3D(cbind(r=r,theta=theta,phi=phi), by=offset$by, ang=offset$ang, sph.in=TRUE, sph.out=sph.out, reverse=reverse)
		#rotate3D(cbind(r=r,theta=theta,phi=phi), by="x", ang=0, sph.in=TRUE, sph.out=sph.out, reverse=reverse)
		}
	else{
		if(!sph.out){
			sph2car(cbind(r=r,theta=theta,phi=phi))
			}
		else{
			cbind(r=r,theta=theta,phi=phi)
			}
		}
	##################################################
	##################################################
	}
