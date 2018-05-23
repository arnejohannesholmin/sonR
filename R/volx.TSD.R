#*********************************************
#*********************************************
#' Calculates the voxel volumes of acoustic instruments implemented in echoIBM (see sonR_implemented()). Usually the volume or area data of one ping is returned, but if muptiple time steps are given in the beam configuration data (specifically for fishery sonars), all time steps are returned if there are changes in the settings between pings.
#'
#' @param data  is the list of inputs variables as returned from read.TSD (including ).
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param var  is a vector of the variables to return. Currently implemented are "volx" for volumes of the voxels and "harx" for horizontal area of the voxels.
#' @param fanWidth  has a number of possible values: (1) "b1": one way beam width. (2) "b2": two way beam width. (3) "fe": beams modeled by rectangular cones with width withing the fan given by the inter-beam angle, and calculated using the equivalent beam angle. This normally causes larged fan width due to overlap between beams.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD labl.TSD NAs strff
#'
#' @export
#' @rdname volx.TSD
#'
volx.TSD<-function(data, esnm="MS70", var=c("volx", "harx"), fanWidth="b2"){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-09-25 - Clean version.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Function for selecting one ping of multiple pings beam configuration data (particularly for fishery sonars):
	getPingBeams = function(x, p){
		x[c("asps","sint","numb","bmmd")] = lapply(x[c("asps","sint","numb","bmmd")], function(xx) if(length(xx)) xx[p] else NULL)
		x[c("lenb","dira","dire","eqba","bwtx","bwty")] = lapply(x[c("lenb","dira","dire","eqba","bwtx","bwty")], function(xx) if(length(xx)>0 && length(dim(xx))==2) xx[,p] else if(length(xx)>0) xx else NULL)
		x$offset$ang = x$offset$ang[p]
		x
		}
	# Extract the name of the acoustic instrument:
	if(!is.null(data$esnm[1])){
		esnm=data$esnm[1]
		}
	
	
	##### Execution #####
	# Get the edges of the voxels, rectangular and/or circular:
	out=list()
	if(sonR_implemented(esnm, "SBE")){
		data = getPingBeams(data, 1)
		data = get.specs.esnm(data, esnm=esnm)
		if(strff("harx",var)){
			#warning("'harx' not available for EK60.")
			out$harx = list()
			}
		if(strff("volx",var)){
			thisdata = volx1D_oneping(data)
			out[names(thisdata)] = thisdata
			}
		}
	else if(sonR_implemented(esnm, "MBE")){
		data = getPingBeams(data, 1)
		data = get.specs.esnm(data, esnm=esnm)
		if(strff("harx",var)){
			#warning("'harx' not available for ME70.")
			out$harx = list()
			}
		if(strff("volx",var)){
			out$volx = NAs(max(data$lenb), max(data$numb))
			# Get the rectangular beams of the fan:
			thisdata = volx2D_oneping(data, fanWidth=fanWidth)
			out$volx[,data$rect] = thisdata$volx
			out$ranges = thisdata$ranges
			out$theta = thisdata$theta
			#out$volx[,data$rect] = volx2D_oneping(data, fanWidth=fanWidth)
			# Get the echosounder beams:
			if(any(!data$rect)){
				out$volx[,!data$rect] = volx1D_oneping(data)$volx
				}
			}
		}
	else if(sonR_implemented(esnm, "MBS")){
		data = getPingBeams(data, 1)
		data = get.specs.esnm(data, esnm=esnm)
		if(strff("harx",var)){
			#warning("'harx' not available for MS70.")
			out$harx = list()
			}
		if(strff("volx",var)){
			thisdata = volx3D_oneping(data)
			out[names(thisdata)] = thisdata
			#out$volx = volx3D_oneping(data)
			}
		}
	else if(sonR_implemented(esnm, "OFS")){
		# If only one time step is given, expand to a matrix for convenience:
		if(length(dim(data$dira))<2){
			presentbeamsvar = intersect(names(data),labl.TSD("rb"))
			data[presentbeamsvar] = lapply(data[presentbeamsvar], function(xx) array(xx, dim=c(length(xx),1)))
			}
		numt = ncol(data$dire)
		# Get the indices of vertical beam mode:
		vertical = data$bmmd==2
		printt(vertical)
		printt(data$dire)
				
		if(strff("harx",var)){
			# Get the pings with identical settings, and consider only the horizontal mode:
			settings = split(which(!vertical), data$dire[1,!vertical])
			
			# If all pings have equal setting, return only one set of voxel horizontal areas:
			if(!any(vertical) && all(data$dire[1]==data$dire[1,])){
				thisdata = getPingBeams(data, 1)
				thisdata = get.specs.esnm(thisdata, esnm=esnm)
				out$harx = harx.SH80_oneping(thisdata)
				}
			else{
				out$harx = NAs(max(data$lenb), max(data$numb), numt)
				for(i in seq_along(settings)){
					# Get the beams data for the first ping of the current settings:
					printt(dim_all(data))
					thisdata = getPingBeams(data, settings[[i]][1])
					printt(dim_all(thisdata))
					thisdata = get.specs.esnm(thisdata, esnm=esnm)
					thisharx = harx.SH80_oneping(thisdata)
					out$harx[seq_len(nrow(thisharx)),,settings[[i]]] = harx.SH80_oneping(thisdata)
					}
				}
			}
		if(strff("volx",var)){
			# Get the pings with identical settings:
			settings = c(split(which(vertical),data$dira[1,vertical]), split(which(!vertical),data$dire[1,!vertical]))
			# If all pings have equal setting, return only one set of voxel volumes:
			if(length(settings)<2){
				thisdata = getPingBeams(data, 1)
				thisdata = get.specs.esnm(thisdata, esnm=esnm)
				thisdata = volx2D_oneping(thisdata, fanWidth=fanWidth)
				out[names(thisdata)] = thisdata
				#out$volx = volx2D_oneping(thisdata, fanWidth=fanWidth)
				}
			else{
				out$volx = NAs(max(data$lenb), max(data$numb), numt)
				for(i in seq_along(settings)){
					# Get the beams data for the first ping of the current settings:
					thisdata = getPingBeams(data, settings[[i]][1])
					thisdata = get.specs.esnm(thisdata, esnm=esnm)
					# Apply again the specification function on the current ping, on order to capture the rotation offset for vertical fans, which are rotated to head horizontally along the vessel direction:
					#thisdata = get.specs.esnm(thisdata)
					thisdata = volx2D_oneping(thisdata, fanWidth=fanWidth)
					out$volx[seq_len(nrow(thisdata$volx)),,settings[[i]]] = thisdata$volx
					out$ranges = thisdata$ranges
					out$theta = thisdata$theta
					#thisdata = volx2D_oneping(thisdata, fanWidth=fanWidth)
					#out$volx[seq_len(nrow(thisdata)),,settings[[i]]] = thisdata
					}
				}
			}
		}
	else{
		warning(paste("Unsupported echo sounder ",esnm," for voxel data",sep=""))
		}
	
	
	##### Output #####
	out
	##################################################
	##################################################
	}
