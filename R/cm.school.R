#*********************************************
#*********************************************
#' Calculates the center of mass of a (simulated) school in 1D, 2D or 3D, or based on acoustic data.
#'
#' @param school  is either a list of x-, y- and z-positions and mass of the fish (named "psxf", "psyf", "pszf" and "mass"), a list of x-, y- and z-positions and volume and acoustic intensity of the voxels (named "psxx", "psyx", "pszx", "volx" and "vbsc" or "mvbs"), or a matrix of columns holding the positions along each dimension.
#' @param mass  is the mass of the fish.
#' @param subset  is a numeric or logical vector defineing which positions to include in the calculation of center of mass.
#' @param na.rm  is used in sum().
#' @param excl.neg  is FALSE if negarive mass values are to be included in the calculation of the centers of mass of the school.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD dim_all NAs zeros
#'
#' @export
#' @rdname cm.school
#'
cm.school<-function(school, mass=1, subset=NULL, na.rm=TRUE, excl.neg=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-02-12 - Clean version.
	# Update: 2011-05-12 - Added the option 'subset', and added the option of calculating the center of mass when voxel positions and volumes, and acoustic intensity data are given.
	# Update: 2011-05-13 - Added the option 'na.rm'.
	# Update: 2013-05-13 - Fixed bug with 'subset'.
	# Last: 2013-08-08 - Fixed bug with school$mass (non-conformable arrays when school$vbsc had multiple time steps).
	########### DESCRIPTION: ###########
	# Calculates the center of mass of a (simulated) school in 1D, 2D or 3D, or based on acoustic data.
	########## DEPENDENCIES: ###########
	# 
	############ VARIABLES: ############
	# ---school--- is either a list of x-, y- and z-positions and mass of the fish (named "psxf", "psyf", "pszf" and "mass"), a list of x-, y- and z-positions and volume and acoustic intensity of the voxels (named "psxx", "psyx", "pszx", "volx" and "vbsc" or "mvbs"), or a matrix of columns holding the positions along each dimension.
	# ---mass--- is the mass of the fish.
	# ---subset--- is a numeric or logical vector defineing which positions to include in the calculation of center of mass.
	# ---na.rm--- is used in sum().
	# ---excl.neg--- is FALSE if negarive mass values are to be included in the calculation of the centers of mass of the school.
	
	
	##################################################
	##################################################
	#collapse = function(school){
	#	dimmass = dim_all(school$mass)
	#	numt = if(length(dimmass)==3) dimmass[3] else 1
	#	len = prod(dimmass)/numt
	#	dim(school$psxf) = c(len, length(school$psxf)/len)
	#	dim(school$psyf) = c(len, length(school$psyf)/len)
	#	dim(school$pszf) = c(len, length(school$pszf)/len)
	#	dim(school$mass) = c(len, numt)
	#	school
	#	}
	
	collapse = function(school){
		dimpos = dim_all(school$psxf)
		numt = if(length(dimpos)==3) dimpos[3] else 1
		len = prod(dimpos)/numt
		dim(school$psxf) = c(len, numt)
		dim(school$psyf) = c(len, numt)
		dim(school$pszf) = c(len, numt)
		dim(school$mass) = c(len, numt)
		school
		}

	
	if(is.list(school)){
		# If segmentation data are present, and 'subset' is empty, consider these as 'subset':
		if(length(subset)==0 && length(school$sgsc)){
			subset=school$sgsc
			}
		# Calculate from voxel positions 'psxx', 'psyx' and 'pszx', and volume backscattering coefficient 'vbsc' (or mean volume backscattering 'mvbs') if present:
		if(!any(length(school$psxx)==0, length(school$psyx)==0, length(school$pszx)==0, length(school$volx)==0, length(school$vbsc)==0)){
			school$psxf=school$psxx
			school$psyf=school$psyx
			school$pszf=school$pszx
			# 2013-08-08: Added c() around school$volx to support multiple time steps in school$vbsc:
			school$mass=school$vbsc*c(school$volx)
			school = collapse(school)
			}
		else if(!any(length(school$psxx)==0, length(school$psyx)==0, length(school$pszx)==0, length(school$volx)==0, length(school$mvbs)==0)){
			school$psxf=school$psxx
			school$psyf=school$psyx
			school$pszf=school$pszx
			# 2013-08-08: Added c() around school$volx to support multiple time steps in school$vbsc:
			school$mass=10^(school$mvbs/10)*c(school$volx)
			school = collapse(school)
			}
		else if(!any( length(school$psxx)==0, length(school$psyx)==0, length(school$pszx)==0)){
			school$psxf=school$psxx
			school$psyf=school$psyx
			school$pszf=school$pszx
			school$mass=zeros(length(school$psxf))+mass
			school = collapse(school)
			}
		# If 'mass' is missing in the list, assign it to be 1 for all positions:
		if(length(school$mass)==0){
			school$mass=zeros(dim_all(school$psxf))+mass
			}
		if(any(school$mass<0,na.rm=TRUE) && excl.neg){
			school$mass[school$mass<0]=0
			warning("Negative masses excluded")
			}
		if(any(school$mass<0,na.rm=TRUE) && !excl.neg){
			warning("Negative mass included. May lead to irregular result")
			}
			
		# Include the logical indexes in 'subset':
		if(length(subset)==0){
			subset=TRUE
			}
		if(is.logical(subset) && length(subset)!=NROW(school$psxf)){
			subset=rep(subset,length.out=NROW(school$psxf))
			}
		
		# Calculate center of mass for all time steps (if arranged along the second dimension):
		if(length(school$psxf)==0  || length(subset)==0){
			NAs(3)
			}
		else{
			drop(cbind(
				colSums(school$psxf[subset,,drop=FALSE]*school$mass[subset,,drop=FALSE],na.rm=na.rm)/colSums(school$mass[subset,,drop=FALSE],na.rm=na.rm),
				colSums(school$psyf[subset,,drop=FALSE]*school$mass[subset,,drop=FALSE],na.rm=na.rm)/colSums(school$mass[subset,,drop=FALSE],na.rm=na.rm),
				colSums(school$pszf[subset,,drop=FALSE]*school$mass[subset,,drop=FALSE],na.rm=na.rm)/colSums(school$mass[subset,,drop=FALSE],na.rm=na.rm)))
			}
		}
	else{
		if(any(mass<0,na.rm=TRUE) && excl.neg){
			mass[mass<0]=0
			warning("Negative masses excluded")
			}
		if(any(mass<0,na.rm=TRUE) && !excl.neg){
			warning("Negative mass included. May lead to irregular result")
			}
		colSums(school*mass,na.rm=na.rm)/sum(mass,na.rm=na.rm)
		}
	##################################################
	##################################################
	}
