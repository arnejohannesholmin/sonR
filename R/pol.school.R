#*********************************************
#*********************************************
#' Cacluates the polarization of a school as one of the following methods: (1) The mean angle deviation between the individual directions and the mean direction of the school, given in radians (Huth and Wissel 1992, Parrish et al 2002). (2) The norm of the average direction unit vectors (Couzin et al 2002). (3) The simple standard deviation of the azimuth and elevation angles of the fish (Holmin et al 2011 (Seattle meeting)).
#'
#' @param x  is either a list of elements 'vlxf', 'vlyf' and 'vlzf', or 'rtzf' and 'rtxf', or a matrix of vectors along the rows.
#' @param type  is a string representing the type of polarization estimate. Currently implemented are: (1) The mean angle deviation between the individual directions and the mean direction of the school (Huth and Wissel 1992, Parrish et al 2002). (2) The norm of the average direction unit vectors (Couzin et al 2002). (3) The simple standard deviation of the azimuth and elevation angles of the fish (Holmin et al 2011 (Seattle meeting)).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD car2sph
#' @importFrom stats sd
#'
#' @export
#' @rdname pol.school
#'
pol.school=function(x,type=c("Huth","Couzin","sd","azimuth","elevation")){
	
	############### LOG: ###############
	# Start: 2011-06-21 - Clean version
	
	########## Preparation ##########
	getHuth=function(x,type="h"){
		meandir=colMeans(x)
		norm_meandir=sqrt(sum(meandir^2))
		if(any(c("a","e") %in% type)){
			# Rotate to have the mean direction in the horizontal plane:
			meandirsph=car2sph(meandir)
			x=rotate3D(x, by="zx", ang=c(-meandirsph[2],pi/2-meandirsph[3]))
			x[,2]=0
			}
		norm_x=sqrt(rowSums(x^2))
		if(any(c("a","e") %in% type)){
			# Recalculate the mean direction:
			meandir=colMeans(x)
			norm_meandir=sqrt(sum(meandir^2))
			}
		# The dot product between two vectors 'a' and 'b' is given by a*b = norm(a) norm(b) cos(theta):
		mean(acos((x %*% meandir) / (norm_meandir * norm_x) ))
		}
	
	# If given as a list, 'x' needs to contain either velocity or angle information:
	if(is.list(x)){
		# Extract angle data:
		x=vl2rt.TSD(x)
		# Combine into matrix:
		if(tolower(substring(type[1],1,1))=="s"){
			x=cbind(x$rtzf,x$rtxf)
			}
		else{
			x=cbind(x$vlxf,x$vlyf,x$vlzf)
			}
		}
	
	
	########## Execution and output ##########
	if(tolower(substring(type[1],1,1))=="h"){
		getHuth(x)
		}
	else if(tolower(substring(type[1],1,1))=="a"){
		warning("The code needs to be checked")
		getHuth(x,"a")
		}
	else if(tolower(substring(type[1],1,1))=="e"){
		warning("The code needs to be checked")
		getHuth(x,"e")
		}
	else if(tolower(substring(type[1],1,1))=="c"){
		# Polarization in [0,1] using the norm of the average direction unit vectors:
		# Get the number of individuals:
		N=nrow(x)
		# Get the length of each individual velocity vector and normalize 'x':
		norm_x=sqrt(rowSums(x^2))
		x=x/norm_x
		# Sum up the direction vectors:
		allx=colSums(x)
		# Take the norm of the resulting velocity and average:
		absallx=sqrt(sum(allx^2))
		absallx/N
		}
	else if(tolower(substring(type[1],1,1))=="s"){
		# Simply return the standard deviation of the azimuth and elevation angle in a list of elements named "sdaz" and "sdel":
		c(sdaz=sd(x[,1]),sdel=sd(x[,2]))	
		}
	else{
		stop("Invalid type (legal are \"Huth\", \"Couzin\" and \"sd\")")
		}
}
