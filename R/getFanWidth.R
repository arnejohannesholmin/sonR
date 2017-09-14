#*********************************************
#*********************************************
#' Calculates the width of a 2-D sonar fan given the specifications in the parameter fanWidth. Only one time step is treated.
#'
#' @param data  is the list of inputs variables as returned from read.TSD (including ).
#' @param fanWidth  has a number of possible values: (1) "b1": one way beam width. (2) "b2": two way beam width. (3) "fe": beams modeled by rectangular cones with width withing the fan given by the inter-beam angle, and calculated using the equivalent beam angle. This normally causes larged fan width due to overlap between beams.
#' @param stretch  is used to stretch the voxels of the ME70 multibeam echosounder in the direction of motion, so that space in between voxels is smoothed out.
#' @param ind  are the indices of the beams to treat (such as data$rect).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom utils tail
#'
#' @export
#' @rdname getFanWidth
#'
getFanWidth <- function(data, fanWidth, stretch=1, ind=NULL, range=0){
	
	# Extract the specified beams:
	numb =length(data$dira)
	ind = ind.expand(data$rect, numb, drop=TRUE)
	beamsVarOfLengthNumb = sapply(data, length) == numb
	dataOriginal = data[beamsVarOfLengthNumb]
	data[beamsVarOfLengthNumb] = lapply(data[beamsVarOfLengthNumb], "[", ind)
	numb <- length(data$dira)
	
	# Rotate the beams to point in the direction of the vessel if vertical is included in the fan. This includes first rotating the coordinate system of the vessel around the z axis so that the x axis coincides with the azimuth angle, and then rotating around the x axis by negative 90 degrees, so point the y axis in the direciton of the fan (downwards):
	vertical = any(data$dire > 170*pi/180)
	if(vertical){
		#dirRotated = rotate3D(cbind(1, data$dira, data$dire), by=paste("z",data$offset$by), ang=c(data$dira[1],data$offset$ang), sph.in=TRUE, sph.out=TRUE, list.out=FALSE)
		dirRotated = rotate3D(cbind(1, data$dira, data$dire), by=data$offset$by, ang=data$offset$ang, sph.in=TRUE, sph.out=TRUE, list.out=FALSE)
		#data$diraOriginal = data$dira
		#data$direOriginal = data$dire
		data$dira = dirRotated[,2]
		data$dire = dirRotated[,3]
		}
	
	#getFanWidth <- function(data, fanWidth, ind=NULL, stretch=1){
	# Get the angles of the rectangular voxels based on the equivalent beam angles and the directions of the beam maxima:
	data$diffdira = abs(diff(data$dira))
	data$diffdira = c(data$diffdira, tail(data$diffdira,1))
	
	#if(length(ind)==0){
	#	ind <- seq_len(N)
	#}
	if(is.numeric(fanWidth)){
		data$diffdire = fanWidth
		}
	else if(strff("ce",fanWidth)){
		warning("Fan width from circular beam using equivalent bema angle not yet implemented")
		}
	else if(strff("ee",fanWidth)){
		warning("Fan width from elliptical beam using equivalent bema angle not yet implemented")
		}
	else if(strff("re",fanWidth)){
		warning("Fan width from rectangular beam using equivalent bema angle not yet implemented")
		}
	else if(strff("fe",fanWidth)){
		# Using the equation for the second angle of a pyramidic cone, from the first angle and the equaivalent beam angle. This equation is found on "http://en.wikipedia.org/wiki/Solid_angle" (Section 2.3 Pyramid, solving for 'a' yealding a = 2*asin( sin(eqba/4) / sin(b/2) )):
		#eqba = 10^(data$eqba[ind]/10)
		eqba = 10^(data$eqba/10)
		arg = sin(eqba/4) / sin(data$diffdira/2)
		# For the verical fan of fishery sonar, the beams closest to vertical orientation are rubbish, and have extremely large equivalent beam angle, resulting in missing across fan angle. These are set to 1 in 'arg' to avoid NAs:
		arg[arg>1] = 1
		data$diffdire = 2*asin( arg )
		}
	else if(strff("b",fanWidth)){
		# Get beam widths as 'bwtx' and 'bwty':
		if(length(data$bwtx)==0){
			data=Simrad_bwt(data)
			}
		# Assume horizontal mode fishery sonar if beam mode is not present:
		if(length(data$bmmd)==0){
			data$bmmd = zeros(numb)
			#warning("The beam mode 'bmmd' not present, set to 0 = horizontal disk")
			}
		# Get the increments in elevation angle:
		data$diffdire = NAs(numb)
		# The beams in the horizontal fan have across-fan beam widths bwty, since the x axis of the coordinate system of the beams is horizontal.
		data$diffdire[data$bmmd==0] = data$bwty[data$bmmd==0]
		# The beams in the vertical fan have across-fan beam widths bwtx, for the same reason.
		data$diffdire[data$bmmd==2] = data$bwtx[data$bmmd==2]
		if(strff("b2",fanWidth)){
			data$diffdire <- data$diffdire / sqrt(2)
			}
		}
	
	# The parameter 'stretch' may be given as a two element vector, where the first element is the stretch value and the second element is the range at which this stretch value is as desired:
	if(length(stretch)==2){
		#stretch=stretch[1]*stretch[2]/r1
		arg=abs(stretch[2])/range * sin(stretch[1]*data$diffdire/2)
		arg[arg>1]=1
		data$diffdire = 2*asin(arg)
		}
	else{
		data$diffdire = data$diffdire * stretch
		}
	
	# Add half the voxel width on each side of the elevation direciton angles:
	data$direExpanded = cbind(data$dire-data$diffdire/2 , data$dire+data$diffdire/2)
	
	# Restore the original beam directions:
	data$dira = dataOriginal$dira
	data$dire = dataOriginal$dire
	
	# Return the data:
	data
	}