#*********************************************
#*********************************************
#' Calculates the midpoints (or edges along the beams, depending on 'rpos') of the voxels of a sonar or an echo sounder, relative to a global coordinate system (usually centered on the sea surface at the first vessel position, with y-axis pointing north and vertical z-axis). Only the heading is considered, as the pitch and roll of the vessel is compensated for by the sonar. If positions are wanted in the coodinate system of the vessel, set 'data' to c(list(psxv=0,psyv=0,pszv=0,rtxv=0,rtyv=0,rtzv=0,psze=0),data).
#' If pitch and/or roll are not compensated for by the sonar, the function has to be altered.
#' The funciton can treat multiple time steps, but not for systems with changing beam configuration between pings, for which the function needs to be run at each separate time step, inputting the beam configuration of one time step at the time (not matrices for variables such as the elevation angle of the direction of the beams 'dire').
#'
#' @param data  is the list of TSD inputs as returned from read.TSD (must contain "psxv", "psyv", "rtzv", "lenb", "numb", "sint", "psze", ("dira","dirl") or ("dire", "dirt"), "asps").
#' @param t  is a single integer giving the time step to treat (in the range [1,number of pings]).
#' @param cs  is the coordinate system of the voxel midpoints or edges of the voxels.
#' @param seabed  is the z-coordinate of the sea bed, usually provided by echo sounder data. Soundbeams reflected from the sea bed or the surface are reflected at the incidence angle.
#' @param rot  is 1 if simple rotation using cosine and sine is to be used, and 2 if the function rotate() is to be used in the rotation. Times for the different methods (tested on MacBook Pro dual 2.8 GHz, 2010-02-09, with other applications running):
#' @param compensation  is a vector of string giving which rotation values that are compensated for in the sonar. Only c("pitch","roll") is available for the current version. Used in soundbeam.TSD.
#' @param ideal  is TRUE to represent the simple case where the speed of sound 'data$asps' is invariant of depth.
#' @param rpos  is a string of two possible values: "midpoint" which specifies that the positions alont the beams are to bo located at midpoints of voxels and "edge" which returns the edges of the voxels including the first edge on the surface of the sonar (length of output is lenb+1 in this case). Abbreviation are allowed.
#' @param drop.out  is TRUE if output should be stripped of dimensions of only one level.
#' @param plot  is TRUE if the plot produced by soundbeam() is to be shown (slows down the calculation).
#' @param ...  allowing for argument passed to tjie function.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom TSD ind.expand zeros numt.TSD
#'
#' @export
#' @rdname soundbeam.TSD
#'
soundbeam.TSD<-function(data, t=1, ind=NULL, cs="g", seabed=-12000, rot=1, compensation=c("pitch","roll"), ideal=TRUE, rpos=c("midp","edge"), drop.out=TRUE, plot=FALSE, ...){
		
	############### LOG: ###############
	# Start: 2009-07-07 - Clean version of the methods ignoring heave.
	# Update: 2009-07-22 - Adding the method accounting for heave.
	# Update: 2009-09-15 - Fixed bug that 'orientations' were interpreted as radians (added the input parameter radians=FALSE) and that heading were treated as counter clockwise (CCW) but should be interpreted as clockwise (CW).
	# Update: 2010-02-09 - Changed method according to the changes in rotate(). Also input chnaged to support list inputs from read.TSD().
	# Update: 2011-09-24 - Added the option 'cs' defineing the coordinate system of the voxel midpoints or edges.
	# Last: 2013-09-12 - Added the option 'ind'.
	
	
	##################################################
	##################################################
	#### Preparation #####
	### Defaults: ###
	# 'data' needs to contain certain elements:
	allpresentideal = ( all(c("rho0","gacc","hpr0","temp","slty") %in% names(data)) | all(c("asps","sint") %in% names(data)) )
	allpresentnotideal = ( all(c("rho0","gacc","hpr0","temp","slty") %in% names(data)) | "isps" %in% names(data) ) & any(c("ihpr","pszc") %in% names(data))
	#allpresent = ( all(c("rho0","gacc","hpr0","temp","slty") %in% names(data)) | "isps" %in% names(data) ) & any(c("ihpr","pszc") %in% names(data))
	
	if(ideal && !allpresentideal){
		if("sint" %in% names(data)){
			data$asps = 1500
			ideal = TRUE
			warning("Average speed of sound defaulted to 1500 m/s, resulting in straight sound beams")
			}
		else{
			warning("The required CTD data are missing, OR \"sint\" missing: All of \"rho0\", \"gacc\", \"hpr0\", \"temp\", and \"slty\" or all of \"asps\" and \"sint\" must be present, in addition to at least one of \"ihpr\" and \"pszc\"")
			}
		}
	if(!ideal && !allpresentnotideal){
		warning("No CTD-data except average speed of sound given, resulting in straight sound beams")
		}
	
	# Check for the presence of vessel dynamic variables when ideal==FALSE:
	cs = cs[1]
	if(!ideal){
		vesselnames = c("psxv","psyv","pszv","rtzv")
		if(!all(vesselnames %in% names(data))){
			stop("'data' must contain the vessel dynamic variables \"psxv\", \"psyv\", \"pszv\" and \"rtzv\" when ideal==FALSE")
			}
		}
	# If ideal==TRUE and cs is "v" set 'pszv' to 0 as a dummy variable used in treatment of 't':
	else if(cs=="v" && !("pszv" %in% names(data))){
		data$pszv = 0
		}
		
	# Length of the inputs:
	l = length(data$pszv)
	# If data$pszv is missing, it is defaulted by zeros:
	if(is.null(data$pszv)){
		data$pszv = zeros(l)
		}
	# Compensation is currently only supported for both pitch (rtzv) and roll (rtxv):
	if(!all(compensation==c("pitch","roll"))){
		stop("Only compensation for both pitch and roll is implemented in this method")
		}
	
	# If the directions of the beams in the spherical coordinatesystem of the vessel are missing, these angles are retrieved from the directions of the beams given in the particular echo sounder sysrem: Currently implemented is data$esnm[1] = "MS70":
	if(any(is.null(data$dira),is.null(data$dire))){
		data = Simrad_dir(data)
		}
	
	# Subsetting the vessel dynamics for conveience:
	data = extractTimeStep(data, t, var=c("dira", "dire", "lenb", "psxv", "psyv", "pszv", "rtzv"))
	#data = extractTimeStep(data, t)
	#data$psxv = data$psxv[t]
	#data$psyv = data$psyv[t]
	#data$pszv = data$pszv[t]
	#data$rtzv = data$rtzv[t]
	
	
	# Subset the beams by 'ind':
	ind = ind.expand(ind, length(data$dira), drop=TRUE)
	data$dira = data$dira[ind]
	data$dire = data$dire[ind]
	numb = length(ind)
	# Dimensional variables:
	lenb = max(data$lenb[ind])
	
	# Create the length of the beams of the output, which will be one longer than the original lengths of the beams if points on the edges are to be returned:
	if(substr(rpos[1],1,1)=="m"){
		lenb_out = lenb
		}
	else if(substr(rpos[1],1,1)=="e"){
		lenb_out = lenb + 1
		}
	else{
		stop("Wrong value of 'rpos', must be one of \"midp\" (representing the midpoints of the voxels along the beams) or \"edge\" (representing the edges of the voxels along beams giving  )")
		}
	
	
	
	# Treatment of 'sonardepth':
	sonardepth = data$psze[1]
	if(sonardepth>0){
		sonardepth = -sonardepth
		warning("'sonardepth' must be non-positive")
		}
	
	# 'uniquedire' is the vetor of the unique elevation angles of the beams. Due to the compensation for pitch and roll perfomed in the echo sounder, beams of the same elevation angle may be replicated, reducing the cpu time of the function to 1/25 for the MS70 sonar:
	uniquedire = unique(data$dire)
	
	
	##### Execution #####
	if(all(is.na(data$rtzv))) {
		warning("All headings (rtzv) are missing, and were set to 0")
		data$rtzv <- 0
	}
	# If ideal==TRUE the soundspeed is assumed to be constant equal to data$asps. In this case simple regularly spaced midpoints are returned:
	if(ideal){
		# Get the ranges:
		r = cbind(soundbeam_range(data, pos=rpos), 0, 0)
		## Radial partitioning:
		#dr = data$asps[1]*data$sint[1]/2
		## Considering the global coordinate system (G), and rotating first around the z-axis by data$dira+data$rtzv, then around the y-axis by data$dire-pi/2, we get a coordinate system with midpoints along the x-axis. Rotating back we get the midpoints in (G):
		#if(substr(tolower(rpos[1]),1,1)=="m"){
		#	r = cbind(c(dr/4,seq_len(lenb_out-1)*dr),0,0)
		#	}
		#else if(substr(tolower(rpos[1]),1,1)=="e"){
		#	r = cbind(c(0,seq_len(lenb_out-1)*dr-dr/2),0,0)
		#	}
		#else{
		#	stop("Wrong value of 'rpos', must be one of \"midp\" (representing the midpoints of the voxels along the beams) or \"edge\" (representing the edges of the voxels along beams giving  )")
		#	}
		
		# The beams are represented in spherical coordinates in the coordinate system of the vessel, which itself is oriented in the global coordinate system. The midpoints/edgepoints of the voxels in the coordinate system of the vessel are obtained by the following steps:
		# 1. Construct a vector of points along the x-axis of the global coordinate system, representing midpoints/edgepoints along the beams.
		# 2. Rotate the global coordinate system by the angles (-dira,-dire) around () 
		
		# z-rotation of the beams will be the azimuth angle with respect to the x-axis:
		theta = data$dira
		# x-rotation of the beams is the elevation angle subtracted pi/2 to equal 0 at the x-axis:
		phi = data$dire-pi/2
		
		if(startsWith(tolower(cs), "v")){
			# Get the positions along the beams:
			xyz = rotate3D(r, ang=cbind(-phi,-theta), by="yz", list.out=TRUE)
			}
		else if(startsWith(tolower(cs), "g")){
			# z-rotation of the beams added the z-rotation of the vessel:
			# c() around data$rtzv avoids the new warning message in R 3.4: "Recycling array of length 1 in vector-array arithmetic is deprecated. Use c() or as.vector() instead."
			theta = data$dira + c(data$rtzv)
			# x-rotation of the beams added the x-rotation of the vessel which is set to 0 according to the compensation for roll and pitch:
			phi = data$dire-pi/2
			# Get the positions along the beams:
			xyz = rotate3D(r, ang=cbind(-phi,-theta), by="yz", list.out=TRUE)
			# Adding the positions of the vessel:
			if(any(length(data$psxv)==0, length(data$psxv)==0)){
				warning("Vessel position data, usually extracted from longitude and latitude information in the NMEA data, are missing. Please add the position data using adds = list(psxv = X-POSITON_VECTOR_OF THE_SAME_LENGTH_AS_TIME_t, psyv = Y-POSITON_VECTOR_OF THE_SAME_LENGTH_AS_TIME_t).")
				numt = numt.TSD(data)
				data$psxv = double(numt)
				data$psyv = double(numt)
				data$pszv = double(numt)
			}
			xyz$x = xyz$x + rep(data$psxv, each=lenb_out*numb)
			xyz$y = xyz$y + rep(data$psyv, each=lenb_out*numb)
			xyz$z = xyz$z + rep(data$pszv + data$psze[1], each=lenb_out*numb)
		}
		else{
			stop(paste("'cs' must be one of \"vessel\" (\"v\") or \"global\" (\"g\"). Was",cs))
			}
		
		# The positions are repeated for all time steps when 'cs' is "v", and when 'cs' is "g" the dimensions are set:
		xyz = lapply(xyz,array,dim=c(lenb_out,numb))
		# Return the ideal positions of the midpoints:
		if(drop.out){
			lapply(xyz,drop)
			return(lapply(xyz,drop))
			}
		else{
			return(xyz)
			}
		}	
	# If ideal==FALSE use soundbeam() to calculate the positions along the beams:
	else{
	
		# The 'x' and 'y' positions of the midpoints of the voxels arranged in [P,I,J,K]-matrices. 
		xyz = list(x=zeros(lenb_out,numb),y=zeros(lenb_out,numb),z=zeros(lenb_out,numb))
		
		thisctd=data[c("rho0","gacc","hpr0","temp","slty","ihpr","pszc","asps")]
		# For each heading of the beams in the center vertical fan (beams number i=13, j=1,...,25) soundbeam() are used to obtaining the midpoints of the 'K' voxels. The sonar is mounted on the port side of the research vessel, which indicates negative x-values:
		for(j in seq_along(uniquedire)){
			# Indexes for and numbers of the beams of equal elevation angle:
			inde = which(data$dire==uniquedire[j])
			linde = length(inde)
			thisx = zeros(lenb_out, linde)
			thisy = zeros(lenb_out, linde)
			thisz = zeros(lenb_out, linde)
			# Convert the elevation angle of the beam to correspond to the definition of 'ang' in soundbeam() (negative for beams heading downwards). That is subtract the elevation angle from pi/2, so that elevation angles larger than pi/2 (below the x-y-plane) give negative 'ang':
			thisang = pi/2-uniquedire[j]
			# Track points along the soundbeams using soundbeam(), and suppress any warnings that may occur:
			thisxyz = suppressWarnings(soundbeam(c(0,data$psze[1]), thisang, ctd=thisctd, seabed=seabed, lenb=lenb, sint=data$sint[1], rpos=rpos, Ncirc=1000, maxlength=750, plot=plot))
			# Two methods are available for the rotation of the soundbeams into each vertical fan for each ping. The first is a simple time saving version that only applies to the situation where the sonar compensates for pitch and roll motion of the vessel, so that only heading is taken into account (here y=0, which leads to the expression from the z-rotation matrix, which also explains the minus sign at y). The second method is a general method using the function rotate3D(), and may be generalized to more than one rotation in the future:
			if(substring(tolower(cs),1,1)=="v"){
				if(rot==1){
					# Setting y=0 in the z-rotation matrix:
					thisx = outer(thisxyz$x,cos(rep(-data$dira[inde])),"*")
					thisy = -outer(thisxyz$x,sin(rep(-data$dira[inde])),"*")
					thisz = rep(thisxyz$y,linde) - data$psze[1]
					}
				else{
					thisxyz = rotate3D(list(x=thisxyz$x,y=double(lenb_out),z=thisxyz$y),by="z",ang=rep(-data$dira[inde]))
					thisx = thisxyz$x
					thisy = thisxyz$y
					thisz = thisxyz$z - data$psze[1]
					}
				}
			else if(substring(tolower(cs),1,1)=="g"){
				if(rot==1){
					thisz = rep(thisxyz$y,linde)
					# Setting y=0 in the z-rotation matrix:
					thisx = outer(thisxyz$x,cos(outer(-data$dira[inde],-data$rtzv,"+")),"*")
					thisy = -outer(thisxyz$x,sin(outer(-data$dira[inde],-data$rtzv,"+")),"*")
					}
				else{
					thisxyz = rotate3D(list(x=thisxyz$x,y=double(lenb_out),z=thisxyz$y), by="z", ang=c(outer(-data$dira[inde],-data$rtzv,"+")))
					thisx = thisxyz$x
					thisy = thisxyz$y
					thisz = thisxyz$z
					}
				# Adding the positions of the vessel:
				thisx = thisx + rep(data$psxv,each=lenb_out*linde)
				thisy = thisy + rep(data$psyv,each=lenb_out*linde)
				thisz = thisz + rep(data$pszv,each=lenb_out*linde)
				}
			else{
				stop(paste("'cs' must be one of \"vessel\" (\"v\") or \"global\" (\"g\"). Was",cs))
				}
			xyz$x[,inde,] = thisx
			xyz$y[,inde,] = thisy
			xyz$z[,inde,] = thisz
			}
		}
		
	
	##### Output #####
	if(drop.out){
		lapply(xyz,drop)
		return(lapply(xyz,drop))
		}
	else{
		return(xyz)
		}
	##################################################
	##################################################
	}
