#*********************************************
#*********************************************
#' This funciton compresses the data stored in 'x' by averaging (either mean, trimmed mean or median) in depth and time bins. This is designed primarily for echosounder beams.
#'
#' @param data  is the list of inputs variables as returned from EKRaw2TSD().
#' @param tres  The time resolution of the compressed data in seconds.
#' @param xres  The sailed distance resolution of the compressed data in meters.
#' @param zres  The depth resolution of the compressed data in meters.
#' @param rres  The range resolution of the compressed data in meters.
#' @param bres  The beam resolution of the compressed data in integer number.
#' @param funvbsc  is the function to apply in the compression, either given as function or as a string, in which case the strings "mean" and "median" represents fast versions of the functions with the corresponding names (sum()/length() and fastMedian(), respectively).
#' @param funt  is the same as funvbsc, but used for averaging vessel data in the new time/distance bins.
#' @param adds  is a list of additional data overriding corresponding variables in 'data'
#' @param split used in psx.TSD().
#' @param skipAngles  is TRUE to discard electircal angles from the data (saves time).
#' @param origin  is either the time index of the origin, or the origin itself, given as c(longitude, latitude).
#' @param z0 is the upper depth of the compression in z direction (vertically), defaulted to 0 (the sea surface).
#' @param drop is TRUE to drop dimensions of the output vbsc (useful if only one frequency is included in the input vbsc (and the dimensions of the vbsc has been dropped)).
#' @param ...  further arguments passed to psx.TSD().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom ccaPP fastMedian
#' @importFrom data.table := data.table key setkeyv
#' @importFrom SimradRaw soundbeam_range
#' @importFrom TSD dim_all global2car labl.TSD NAs strff utim.TSD utim2mtim zeros
#' @importFrom utils tail head
#'
#' @export
#' @rdname compr.TSD
#'
compr.TSD <- function(data=NULL, tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("median","mean"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, z0=0, drop=FALSE, ...){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2011-09-25 - Clean version.

	
	##################################################
	##################################################
	# Function stolen from http://stackoverflow.com/questions/7735647/replacing-nas-with-latest-non-na-value. Repeats the last non-missing value. Used when filling in data in empty time intervals:
	repeatNonNA <- function(x){
		olddim = dim(x)
		ind = which(!is.na(x))
		# If it begins with a missing, add the first position to the indices:
		if(is.na(x[1])){
			ind = c(1, ind)
		}
		# Repeat the values at these indices diffing the indices + length yields how often :
		x = rep(x[ind], times=diff(c(ind, length(x) + 1) ))
		dim(x) = olddim
		x
	}		
	# Function used for applying the function 'funt':
	applyFunt <- function(x, tindex, funt, numt=NULL){
		applyFuntOne = function(xx){
			# Expand the indices to the array:
			olddim = dim_all(xx)
			len = length(xx)
			n1 = len/length(tindex)
			ind1 = rep(tindex, each=n1)
			ind2 = rep(seq_len(n1), length(tindex))
			
			# Apply the function 'funt' to the time intervals:
			out = c(tapply(xx, list(ind2, ind1), function(xxx) if(is.character(xxx)) head(xxx,1) else funt(xxx)))
			utindex = unique(tindex)
			dim(out) = c(olddim[-length(olddim)], length(utindex))
			
			# If numt is given larger than the number of unique time steps, return an expanded array with NAs at missing time intervals:
			if(length(numt) && numt>length(utindex)){
				outWithNA = NAs(olddim[-length(olddim)], numt)
				# Special care for one-dimensional vectors:
				if(length(olddim)==1){
					outWithNA[utindex] = out
				}
				else{
					outWithNA[, utindex] = out
				}
				out = repeatNonNA(outWithNA)
			}
			out
		}
		lapply(x, applyFuntOne)
	}
	
	if(length(data$vbsc)==0){
		return()
	}
	if(any(sapply(c(tres, xres, zres, rres, bres), length)>0)){
		# Merge the 'adds' and the data:
		if(length(adds)>0){
			data[names(adds)] = adds
		}
		data$utim = utim.TSD(data[labl.TSD("t")])

		# Get the dimensions of the vbsc
		dimvbsc = c(max(data$lenb), data$numb[1], length(data$utim))
		outdim = dimvbsc
		rindex = NULL
		tindex = NULL
		bindex = NULL
		origin = c(data$lonv[origin], data$latv[origin])
		if(is.null(data$psxv) || is.null(data$psyv)){
			data$psyv = global2car(cbind(c(data$lonv), c(data$latv)), origin=origin)
			data$psxv = data$psyv[, 1]
			data$psyv = data$psyv[, 2]
		}
		
		# Get the names of each category:
		beamsnames = setdiff(intersect(names(data), labl.TSD("b")), labl.TSD("t"))
		vesselnames = intersect(names(data), labl.TSD("v"))
		timenames = intersect(names(data), labl.TSD("t"))
		
		if(is.character(funvbsc)){
			if(strff("mean", funvbsc[1])){
				funvbsc = function(xx) sum(xx, na.rm=TRUE)/length(xx)
			}
			else if(strff("median", funvbsc[1])){
				funvbsc = fastMedian
			}
			else{
				funvbsc = get(funvbsc[1])
			}
		}
		
		if(is.character(funt)){
			if(strff("mean", funt[1])){
				funt = function(xx) sum(xx, na.rm=TRUE)/length(xx)
			}
			else if(strff("median", funt[1])){
				funt = fastMedian
			}
			else{
				funt = get(funt[1])
			}
		}
		
		
		###########################################
		########## 1. Compress in range: ##########
		###########################################
		# Compression in range 
		if(length(rres)>0){
			rindex = rep(soundbeam_range(data[c("lenb", "rres", "asps", "sint")], "mid"), prod(tail(dimvbsc,2)))
			# Split the data into range bins:
			if(length(rres)==1){
				if(length(data$rres)==0){
					# Get range resolution:
					data$rres = data$asps[1] * data$sint[1]/2
				}
				# Change rres and lenb to generate range intervals:
				rindexintervals = soundbeam_range(data[c("lenb", "rres", "asps", "sint")], pos="grid", adds=list(lenb=rep(ceiling(data$lenb[1] * data$rres[1] / rres), data$numb[1]), rres=rep(rres, l=data$numb[1])))
				#rangeintervals = seq(floor(min(rindex)/rres), ceiling(max(rindex)/rres)) * rres
			}
			else{
				rindexintervals = rres
			}
			rindex = findInterval(rindex, rindexintervals, all.inside=TRUE)
		}
		### 2. ... or compress in depth: ###
		else if(length(zres)>0){
			rindex = psx.TSD(data, pad=TRUE, split=split, ...)$psz
			# Split the data into depth bins:
			if(length(zres)==1){
				rindexintervals = seq(floor(min(rindex, na.rm=TRUE)/zres), z0) * zres
			}
			else{
				rindexintervals = zres
			}
			rindex = findInterval(rindex, rindexintervals, all.inside=TRUE)
			# Reverse the intervals, since depth is negative:
			rindex = length(rindexintervals) - rindex
			# Modify the transducer depth to z0:
			if(length(dim(data$psze))==2){
				data$psze = matrix(z0, nrow=1, ncol=dimvbsc[3])
			}
			else{
				data$psze = rep(z0, length.out=dimvbsc[3])
			}
			
		}
		# Alter the dimensions, and change the length of beams and range resolution:
		if(length(rindex)){
			outdim[1] = length(rindexintervals)-1
			# Change the length of the beams and the radial resolution
			data$lenb = matrix(length(rindexintervals)-1, nrow=dimvbsc[2], ncol=dimvbsc[3])
			data$rres = rep(tail(diff(rindexintervals),1), l=dimvbsc[3])
		}
		else{
			rindex = rep(seq_len(dimvbsc[1]), prod(dimvbsc[-1]))
		}
		
		##########################################
		########## 3. Compress in time: ##########
		##########################################
		if(length(tres)>0){
			vec = data$utim
			# Split the data into time bins:
			if(length(tres)==1){
				tindexintervals = seq(floor(min(vec, na.rm=TRUE)/tres), ceiling(max(vec, na.rm=TRUE)/tres)) * tres
			}
			tindex = findInterval(vec, tindexintervals, all.inside=TRUE)
		}
		### 4. ... or compress over sailed distance: ###
		else if(length(xres)>0){
			# Convert to meters:
			#tindex = data$utim[1] + (data$sadv - min(data$sadv)) * 1852
			# The above added time to distance, which does not make sense:
			vec = data$sadv * 1852
			# Split the data into sailed distance bins:
			if(length(xres)==1){
				tindexintervals = seq(floor(min(vec, na.rm=TRUE)/xres), ceiling(max(vec, na.rm=TRUE)/xres)) * xres
			}
			tindex = findInterval(vec, tindexintervals, all.inside=TRUE)
		}
		# Alter the dimensions, and subset the beams and vessel data:
		if(length(tindex)){
			
			# Add log, lon, lat and time start and end:
			start_t <- which(!duplicated(tindex))
			end_t <- c(start_t[-1], length(tindex))
			
			#data$start_t <- start_t
			#data$end_t <- end_t
			
			data$log1 <- data$sadv[start_t]
			data$lon1 <- data$lonv[start_t]
			data$lat1 <- data$latv[start_t]
			data$utm1 <- data$utim[start_t]
			
			data$log2 <- data$sadv[end_t]
			data$lon2 <- data$lonv[end_t]
			data$lat2 <- data$latv[end_t]
			data$utm2 <- data$utim[end_t]
			
			
			# Get the unique time steps indices:
			outdim[3] = length(tindexintervals) - 1
			
			# Extract beams data. Here the use of numt either fills the empty time intervals with NAs:
			if(length(dim(data$freq))==2){
				beamsToBeChanged = sapply(data[beamsnames], function(xx) if(length(xx)>0) tail(dim_all(xx),1)==dimvbsc[3] else FALSE)
				data[beamsnames][beamsToBeChanged] = applyFunt(data[beamsnames][beamsToBeChanged], tindex, funt, numt=outdim[3])
			}
			# Extract vessel data:
			thesevesselnames = setdiff(vesselnames, timenames)
			data[thesevesselnames] = applyFunt(data[thesevesselnames], tindex, funt, numt=outdim[3])
			#data[thesevesselnames] = lapply(data[thesevesselnames], function(xx) tapply(xx, tindex, funt))
			
							# Extract the time information using all intervals:
							diffTindexintervals = diff(tindexintervals)
							temp <- tindexintervals[-length(tindexintervals)] + diffTindexintervals/2
							if(length(tres)>0){
								data$utim <- temp
								#data$sadv <- 
							}
							### 4. ... or compress over sailed distance: ###
							else if(length(xres)>0){
								data$utim = (data$utm1 + data$utm1) / 2
								data$sadv <- temp
							}
							
							# Get the number of time steps in each compressed time bin:
							data$nmtc = zeros(outdim[3])
							data$nmtc[unique(tindex)] = table(tindex)
			
			#data$utim = tindexintervals[presentTindexintervals] + diffTindexintervals[presentTindexintervals]/2
			##data$utim = tindexintervals[-length(tindexintervals)] + diff(tindexintervals)/2
			
			# Add Matalb time (but remove this first and generate from the utim):
			data$mtim <- NULL
			data$mtim = mtim.TSD(data)
			data$indt = NULL
			#data$tnxi = tindexintervals
			#data$tndx = tindex
			tindex = rep(tindex, each=prod(head(dimvbsc,2)))
		}
		else{
			tindex = rep(seq_len(dimvbsc[3]), each=prod(dimvbsc[1:2]))
		}
		
		#############################################
		########## 5. Compress over beams: ##########
		#############################################
		oldindi = seq_len(data$numb[1])
		# We need this to separate between beams:
		bindex = rep(rep(oldindi, each=dimvbsc[1]), dimvbsc[3])
		if(length(bres)>0){
			# Split the data into beam bins:
			if(length(bres)==1){
				bindexintervals = 0.5 + seq(0, ceiling(max(oldindi)/bres)) * bres
			}
			else{
				bindexintervals = bres
			}
			# Set the new dira and dire:
			olddiredirainnew = findInterval(oldindi, bindexintervals, all.inside=TRUE)
			
			# Treat beams data given independent of time (equal for all pings):
			if(length(dim(data$dira))==0){
				beamsToBeChanged = sapply(data[beamsnames], function(xx) length(xx)==max(sapply(data[beamsnames], length)))
				data[beamsnames][beamsToBeChanged] = lapply(data[beamsnames][beamsToBeChanged], function(xx) by(xx, olddiredirainnew, funt))
				data$numb = length(data$dira)
				bindex = findInterval(bindex, bindexintervals, all.inside=TRUE)
				outdim[2] = length(bindexintervals)-1
			}
			# Treat beams data given per ping:
			else{
				beamsToBeChanged = sapply(data[beamsnames], function(xx) NROW(xx)==max(sapply(data[beamsnames], NROW)))
				data[beamsnames][beamsToBeChanged] = lapply(data[beamsnames][beamsToBeChanged], function(xx) do.call(rbind, by(xx, olddiredirainnew, function(yy) apply(yy, 2, funt))))
				data$numb = apply(data$dira, 2, length)
				dim(data$numb) <- c(1, length(data$numb))
				bindex = findInterval(bindex, bindexintervals, all.inside=TRUE)
				outdim[2] = length(bindexintervals)-1
			}
		}
		
		
		# Apply the compression
		if(!skipAngles && length(data$angt)){
			DT = data.table(vbsc=c(data$vbsc), angt=c(data$angt), angl=c(data$angl))
		}
		else{
			data$angt = NULL
			data$angl = NULL
			DT = data.table(vbsc=c(data$vbsc))
		}
		
		# Add the indices:
		presentKeys = NULL
		if(length(tindex)){
			DT[,tindex:=tindex]
			presentKeys = c(presentKeys, "tindex")
		}
		if(length(bindex)){
			DT[,bindex:=bindex]
			presentKeys = c(presentKeys, "bindex")
		}
		if(length(rindex)){
			DT[,rindex:=rindex]
			presentKeys = c(presentKeys, "rindex")
		}
		
		
		setkeyv(DT, presentKeys)
		
		data$vbsc = NAs(outdim)
		temp = DT[, funvbsc(vbsc), by=key(DT)]
		temp2 = as.matrix(temp[, rev(presentKeys), with=FALSE])
		data$vbsc[temp2] = temp$V1
		
		
		if(!skipAngles && length(data$angt)){
			warning("Electric angles are deprecated")
			
			data$angt = NAs(outdim)
			temp = DT[, funvbsc(angt), by=key(DT)]
			data$angt[as.matrix(temp[, rev(presentKeys), with=FALSE])] = temp$V1
			
			data$angl = NAs(outdim)
			temp = DT[, funvbsc(angl), by=key(DT)]
			data$angl[as.matrix(temp[, rev(presentKeys), with=FALSE])] = temp$V1
		}
		
		data[c("psxx", "psyx", "pszx")] = psx.TSD(data, pad=TRUE, split=split, ...)
	}
	else{
		warning("No compression")
	}
	
	if(drop){
		data$vbsc <- drop(data$vbsc)
	}
		
	
	##### Output #####
	invisible(data)
	##################################################
	##################################################
	}
