#*********************************************
#*********************************************
#' Computes the path of a ray of sound for the given temperature, salinity and depth/pressure values. The sound beam is calculated in the x-y-plane, where y corresponds to z in three dimensions. If unequal in length, the shorter ones of the inputs are recycled. A list of the following elements is returned Used in soundbeam.TSD().
#'
#' @param start  is the transducer position.
#' @param ang  is the initial angle of the sound beam, defined on the unit circle in the interval (-pi/2,pi/2) (for downward starting angle ang is negative).
#' @param ctd  is the conductivity-temperature-depth data given as a list of suitable names.
#' @param seabed  is the z-position of the sea bed.
#' @param time  is the vector of time points at which the positions along the sound beam are recorded. If not given, 'time' is calculated from 'lenb' and 'sint'.
#' @param w  is the vector of positions along the sound beam, overriding 'time' or 'lenb' and 'sint'.
#' @param lenb  is the number of partitions radially of equal distance 'sint' in time. 
#' @param sint  is the time length of the sound pulse. 'lenb' and 'sint' are used specifically when calculating midpoints of voxels of sonars, where the sound beam travels out to an object and is reflected back to the sonar. Thus the effective time used is dt=sint/2 (see the master thesis of Arne Johannes Holmin page 30-33). The first time point starts at dt/4, and the remaining time points are dt*k, k=2,3,...,lenb. (See the master thesis of Arne Johannes Holmin page 31-32.)
#' @param rpos  is a string of two possible values: "midpoint" which specifies that the positions alont the beams are to bo located at midpoints of voxels and "edge" which returns the edges of the voxels including the first edge on the surface of the sonar (length of output is lenb+1 in this case). Abbreviation are allowed.
#' @param Ncirc  are the maximum number of circles from which the soundbeams are reconstructed.
#' @param maxlength  is the maximum length of the beam, restricting the number of calculations of circles.
#' @param plot  is TRUE if the sound beams should be plotted (time consuming).
#' @param x  is a vector of the x-positions of the sound beam.
#' @param y  is a vector of the y-positions of the sound beam.
#' @param l  is a vector of the tracked positions along the sound beam.
#' @param r  is a vector of the origins of the circles along which the sound beam is traced for each sound speed layer.
#' @param ang  is a matrix of two columns representing the vectors of the start and end angles of the circle segments along which the sound beam is traced for each sound speed layer, as defined on the unit circle.
#' @param origin  is a matrix of two columns representing the x- and y-positions of the points linking the circle segments.
#' @param pos  is a vector of the x-positions of the sound beam.
#' @param extrapolated  is a logical vector which is TRUE for all circle segments located in sound speed layers extrapolated from the ctd-data, in the case that the ctd measurements do not extend the entire water column.
#' @param totlength  is the total length of the calculated path along which the sound beam positions are tracked.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD zeros
#' @importFrom graphics lines points
#'
#' @export
#' @rdname soundbeam
#'
soundbeam<-function(start=c(0,0), ang, ctd, seabed=-12000, time=NULL, w=NULL, lenb=1500, sint=5.12e-4, rpos=c("midpoint","edge"), Ncirc=1000, maxlength=Inf, plot=TRUE, Pain=TRUE){
	
	############### LOG: ###############
	# Start: 2008-02-23 - Finished.
	# Update: 2009-06-12 - Changed to support list and matrix input. Fundamental change in method from the old soundbeam(), which used the old function propagate() for the tracking along the path of the sound beam (see the "-unused" directory for old version).
	# Update: 2009-07-01 - Changed to suppert tracking along the sound beams by equally separated time intervals as well as length intervals. 
	# Update: 2009-07-12 - Changed in method when adjusting the depth 'D' and speed 'speed' to match the range c(seabed,0), by using the function extrapolate.matrix().
	# Update: 2009-07-23 - Fixed bug occuring when the start point is located in the interior of a layer, and the first circle segment does not cross into the next layer. The function structure is changed to defining all variables at the start of the "execution" section, and initializing all variables before the while loop, and then updating the variables at the end of the while loop. Description updated and systemized.	
	# Last: 2009-09-04 - Fixed bug for the special case when ang=0 and the beam starts in the border between two layers.	
	
	##################################################
	##################################################
	##### Preparation #####
	# 'ctd' needs to contain certain elements:
	allpresent=( all(c("rho0","gacc","hpr0","temp","slty") %in% names(ctd)) | "isps" %in% names(ctd) ) & any(c("ihpr","pszc") %in% names(ctd))
	if(!allpresent){
		if(!("asps" %in% names(ctd))){
			stop("'ctd' need to contain the element \"isps\" (instantaneous speed of sound) or all of \"rho0\" (mass density), \"gacc\" (gravitational accelleration), \"hpr0\" (hydrostatic pressure at sea level), \"temp\" (teperature), \"slty\" (salinity), in addition to one of \"ihpr\" (instantaneous hydrostatic pressure), \"pszc\" (depth in negative meters). If not fulfilled, \"asps\" (average speed of sound) may be given, resulting in straight sound beams.")
			}
		else{
			warning("Constant speed of sound ('asps') used when calculating the sound beams")
			}
		}
	
	# Internal function finding the radii of the circle segments:
	radius<-function(ang,speed,eta){
		speed/eta/sin(ang)
		}
	
	# Sea bed must be below sea level:
	if(seabed>=0){
		warning("Non-negative 'seabed'. Interpreted as negative")
		seabed=-seabed
		}
	if(allpresent){
		speed=speedofsound(ctd,Pain=Pain)
		D=getzfromctd(ctd,Pain=Pain)
		}
	else{
		D=c(1,0,start[2]-sqrt(.Machine$double.eps),seabed,seabed-1)
		speed=c(Inf,ctd$asps,ctd$asps,ctd$asps,Inf)
		}
	
	# Stripping 'ctd' of duplicate rows:
	D_speed=cbind(D,speed)
	D_speed=unique(D_speed)
	# Order 'D' and 'speed':
	orderD=order(D_speed[,1])
	rangeD=D_speed[c(1,length(D_speed[,1])),1]
	# Extrapolate D_speed to the range c(seabed,0):
	D_speed=extrapolate.matrix(D_speed[orderD,,drop=FALSE],c(seabed,0),along=1)
	# Adding Inf-values below the sea bed and above the sea surface, causing the sound beam to reflect at these points:
	D_speed=rbind(c(seabed-1,Inf),D_speed,c(1,Inf))
	
	# Splitting up into depth 'D' and speed 'speed':
	D=D_speed[,1]
	speed=D_speed[,2]
	
	# Check if any of the water layers specified by the ctd-data have speed of sound invariant of depth:
	diffspeed=diff(speed)
	# Add or subtract .Machine$double.eps^(1/5) from the consecutive sound speed values that are equal, so that infinite radii is avoided:
	if(any(diffspeed==0,na.rm=TRUE)){
		speed[diffspeed==0]=speed[diffspeed==0]+.Machine$double.eps^(1/5)*rep(c(-1,1),l=sum(diffspeed==0))
		}
	
	# 'extrapolated' records which speed values that are extrapolated to the sea surface and to the sea bed:
	extrapolated=D>rangeD[1] | D<rangeD[2]
		
	# 'nlayers', 'diffD', 'diffspeed' and 'midspeed' defines the speed properties of the layers, and 'eta' is the speed gradient of the layers
	nlayers=length(D)
	diffD=diff(D)
	diffspeed=diff(speed)
	midspeed=speed[-1]-diffspeed/2
	eta=diffspeed/diffD
	# 'midspeed' needs to be adjusted to the closest value at the ends, to ensure that radius() works (since the speed has been given Inf values at the ends):
	midspeed[c(1,nlayers-1)]=midspeed[c(2,nlayers-2)]
	
	# The start point cannot be located above the surface:
	if(start[2]>0){
		start[2]=-start[2]
		}
	
	# If vertical starting angle is specified, add a tiny value to the starting angle (cheating, but simple):	
	if(ang %% (pi/2) == 0){
		ang=ang+sqrt(.Machine$double.eps)
		}	
	
	# If the value of 'ang' is given outside of the interval [-pi/2,pi/2], it is mapped on to this interval by using asin(sin(ang)):
	mirrored=FALSE
	if(ang<(-pi/2) || pi/2<ang){
		ang=asin(sin(ang))
		mirrored=TRUE
		}
	
	
	##### Execution #####
	
	### The variables used in the calculation: ###
	# (1,2,3,4,5) are arrays that are updated at each step and returned at the end (A):
	# (6,7,8,9,10,11,12,13,14) are dummy variables, only used in the calculation (D):
	
	# (1) origin (A)
	# (2) r (A)
	# (3) angles (A)
	# (4) linkpoints (A)
	# (5) layer (A)
	# (6) etaneg (D)
	# (7) lh (D)
	# (8) fromlast (D)
	# (9) h (D)
	# (10) angdev (D)
	# (11) up (D)
	# (12) quadrant2or3 (D)
	# (13) notcross (D)
	# (14) totlength (D)
	
	# 'origin' are the origins of the circle segments:
	origin=array(0,dim=c(Ncirc+1,2)) # (1)
	# 'r' are the radii of the circle segments:
	r=double(Ncirc+1) # (2)
	# 'angles' are the angles of the circle segmetns, as defined on the unit circle:
	angles=array(0,dim=c(Ncirc+1,2)) # (3)
	# 'linkpoints' are the start points and end points of the circle segments (which may be found by the 'origin', 'r' and 'angles'):
	linkpoints=array(0,dim=c(Ncirc+1,2)) # (4)
	# 'layer' are the numbers of the layers in which the sound beam is located:
	layer=double(Ncirc+1) # (5)
	# 'etaneg' is TRUE if 'eta' is negative in the current layer. If so, pi must be added to the angles in 'angles', because this imply one of the quadrants 3 or 4:
	etaneg=FALSE # (6)
	# 'lh' is the distance vector from the 'origin' to the current start point of the circle segment:
	lh=double(2) # (7)
	# 'fromlast' is the distance from the start point z-value to the border which the sound beam appears to have crossed. Only used at the initial step to assure that the value 'h' (se next variable) is compensated for the position of the sonar when not on exactly one of the depth values (!any(start[2])==D)):
	fromlast=0 # (8)
	# 'h' is the height of the current circle segment, i.e. the range in z-value. If the start point is located in the interior of a layer, the absolute value of 'fromlast' is added:
	h=0 # (9)
	# 'angdev' is the difference between the current start angle and vertical, used if the circle segment never crosses into the next layer:
	angdev=0 # (10)
	# 'up' is TRUE if the sound beam points upwards:
	up=FALSE # (11)
	# 'quadrant2or3' is TRUE if the start of the circle segment is located in the 2. or the 3. quadrant:
	quadrant2or3=FALSE # (12)
	# 'notcross' is TRUE if the current circle segment does not cross into the next layer, implying that the circle segment is shaped as an arc:
	notcross=FALSE # (13)
	# 'totlength' is the total length of the sound beam at the current step:
	totlength=0 # (14)
	
	### Initial values: ###
	# If ang=0 and the beam starts in the border between two layers, a tiny distance is subtracted from start[2] (2009-09-04):
	if(any(start[2]==D) && ang%%pi==0){
		smalladd=-sqrt(.Machine$double.eps)
		start[2]=start[2]+smalladd
		}
	linkpoints[1,]=start # (4)
	layer[1]=sum(D<start[2]) # (5)
	etaneg=eta[layer[1]]<0 # (6)
	# The first value of 'angles' is added pi/2 to transform from direction of sound beam to angle of the circle segment:
	angles[1,1]=ang+pi/2+etaneg*pi # (3)
	# Utilizing the internally defined radius() function:
	r[1]=radius(angles[1,1],midspeed[layer[1]],eta[layer[1]]) # (2)
	
	# Obtaining the first point of 'origin':
	lh=r[1]*c(cos(angles[1,1]),sin(angles[1,1])) # (7)
	origin[1,]=linkpoints[1,]-lh # (1)
	# Obtaining the heigth of the circle segment, to check if it crosses into the next layer or simply moves back into the previous layer:
	fromlast=linkpoints[1,2]-D[layer[1]+etaneg] # (8)
	h=r[1]*(1-abs(sin(angles[1,1])))+abs(fromlast) # (9)
	# 'angdev' is used when obtaining the next 'angles' value, if the circle segment does not cross into the next layer:
	angdev=angles[1,1]%%pi-pi/2 # (10)
	# It can be showh that 'angdev' is positive if the sound beam moves upwards:
	up=angdev>0 # (11)
	# The sound beam is in the 2. or 3. quadrant if it moves upwards and eta in the current layer is positive, implying downwards refraction, or if it point downwards and eta in the current layer is negative, implying upwards refraction:
	quadrant2or3=(etaneg!=up) # (12)
	# The sound beam does not cross into the next layer if the heigth 'h' of the circle segment is smaller than the height of the current layer and if the sound beam is in the 2. or 3. quadrant:
	notcross = h<=diffD[layer[1]] & quadrant2or3 # (13)
	# Adjusting 'angdev' and 'lh' if the first circle does not cross into the next layer (notcross==TRUE) and the start postition is in the interior of a layer (!any(start[2])==D)).
	if(notcross){
		
		### Error in if (notcross) { : missing value where TRUE/FALSE needed
		
		# 'angdevright' is the angle span of the part of the circle segment to the right of the vertical radius line:
		angdevrigth=pi/2-asin((D[layer[1]+etaneg]-origin[1,2])/r[1])%%(2*pi)%%pi
		# Adding the x-propagation of the circle path to the right of the vertical radius line (the minus sign and the abs() term is to assure that a negative value is added, as lh[1] is negative when notcross==TRUE):
		lh[1]=(lh[1]-r[1]*sin(abs(angdevrigth)))/2 # (7)
		# Adding the angle of the circle path to the right of the vertical radius line:
		angdev=(angdevrigth+angdev)/2 # (10)
		}
		
	### While loop constructing the path of the sound beam: ###
	i=1
	while(1){
		
		# If the the sound beam does not cross into the next layer, the direction up/down changes:
		if(notcross){
			# Moving back to the previous layer:
			layer[i+1]=layer[i]+(-1)^up # (5)
			linkpoints[i+1,2]=D[layer[i]+!up] # (4)
			# Simply adding the corde:
			linkpoints[i+1,1]=linkpoints[i,1]-2*lh[1] # (4)
			# Subtracting 2 times the difference between the angle at the start of the current circle segment and vertical:
			angles[i,2]=angles[i,1]-2*angdev # (3)
			}
		else{
			# Moving to the next layer:
			layer[i+1]=layer[i]+(-1)^!up # (5)
			linkpoints[i+1,2]=D[layer[i]+up] # (4)
			# Using the information in 'linkpoints[i+1,2]' to find the angles of the current circle segment. The use of 'quadrant2or3' is needed to compensate for angles being placed in the wrong quadrant by asin():
			angles[i,2]=(-1)^quadrant2or3 * asin((linkpoints[i+1,2]-origin[i,2])/r[i])+pi*quadrant2or3 # (3)
			angles[i,2]=angles[i,2]%%(2*pi) # (3)
			linkpoints[i+1,1]=origin[i,1]+r[i]*cos(angles[i,2]) # (4)
			}
		
		# If the maximum number of circle segments is reached, or if the maximum length of the sound beam is exceeded, the fucntions stops and returns the results:
		if(i>=Ncirc || totlength>maxlength){
			break
			}
		# Updating 'totlength':
		totlength=totlength+r[i]*abs(c(angles[i,1]-angles[i,2])) # (14)
		
		### Updating variables: ###
		i=i+1
		# Updating 'etaneg', 'angles', 'r', 'lh', 'origin', 'h', 'angdev', 'up', 'quadrant2or3' and 'norcross':
		etaneg=eta[layer[i]]<0 # (6)
		angles[i,1]=etaneg*pi+angles[i-1,2]%%pi # (3)
		r[i]=radius(angles[i,1],midspeed[layer[i]],eta[layer[i]]) # (2)
		lh=r[i]*c(cos(angles[i,1]),sin(angles[i,1])) # (7)
		origin[i,]=linkpoints[i,]-lh # (1)
		h=r[i]*(1-abs(sin(angles[i,1]))) # (9)
		angdev=angles[i,1]%%pi-pi/2 # (10)
		up=angdev>0 # (11)
		quadrant2or3=(etaneg!=up) # (12)
		notcross = h<=diffD[layer[i]] & quadrant2or3 # (13)
		}
	
	
	##### Return #####
	# Stripping the outputs of zeros, if 'maxlength' kicked in before Ncirc:
	valid=1:(i-1)
	origin=origin[valid,,drop=FALSE] # (1)
	r=r[valid] # (2)
	angles=angles[valid,,drop=FALSE] # (3)
	linkpoints=linkpoints[c(valid,i),,drop=FALSE] # (4)
	layer=layer[valid] # (5)
	# Locating points of equal distance along the sound beam. The first midpoint is on the sonar surface while the rest are distributed in sint/2 intervals, but since it is unwanted to have the first medpoint on the sonar surface, we position it in the middle of the first voxel (0,sint/4) -> sint/8 = 0.125*sint:
	warning("soundbeam.range() should be implemented in soundbeam()")
	if(is.null(time)){
		if(substr(rpos[1],1,1)=="m"){
			dt=c(0.125,0.375,rep(0.5,length.out=lenb-2))*sint
			}
		else if(substr(rpos[1],1,1)=="e"){
			dt=c(0,0.25,rep(0.5,length.out=lenb-1))*sint
			}
		else{
			stop("Wrong value of 'rpos', must be one of \"midp\" (representing the midpoints of the voxels along the beams) or \"edge\" (representing the edges of the voxels along beams giving  )")
			}
		}
	else{
		dt=c(time[1],diff(time))
		}
	# If mirrored==TRUE the points returned by track.circle() are mirrored around the line x=start[1], in which case the plotting is taken outside of track.circle():
	if(mirrored){
		xy=track.circle(origin=origin,r=r,ang=angles,dw=dt,plot=FALSE,speed=midspeed[layer])
		# Mirror the x-value2s:
		if(mirrored){
			xy$x=start[1]-(xy$x-start[1])
			}
		# Plot the points and circles (taken from track.circle()):
		if(plot){
			dens=100
			lx=length(r)
			circles=zeros(lx,dens,2)
			for(i in 1:lx){
				circles[i,,]=circle(c(origin[i,1],origin[i,2]),r[i],seq(angles[i,1],angles[i,2],length.out=dens))
				}
			circles[,,1]=start[1]-(circles[,,1]-start[1])
			plot(NULL,xlim=range(xy$x,circles[,,1]),ylim=range(xy$y,circles[,,2]))
			for(i in 1:lx){
				lines(circles[i,,])
				}
			points(xy$x,xy$y,col=2,pch="*",cex=1)
			}
		}
	else{
		xy=track.circle(origin=origin,r=r,ang=angles,dw=dt,plot=plot,speed=midspeed[layer])
		}
	# Output:
	list(x=xy$x,y=xy$y,l=xy$l,r=r,ang=angles,origin=origin,pos=linkpoints,extrapolated=extrapolated[layer],totlength=totlength)
	##################################################
	##################################################
	}
