#*********************************************
#*********************************************
#' Merges time steps of acoustic data located in the neighbor directory "temp" of the "tsd"-directory given by 'event'.
#'
#' @param x  a list containing the possibly compressed acoustic data (vbsc, vxIX).
#' @param pad  is TRUE to pad with missing values to make the list elements fit in the dimension up to the last dimension.
#' @param split  Logical: if TRUE, split each list element to have last dimenstion equal to 1.
#' @param filename  is the filename of the file containing the data, used in warning messages.
#' @param t		is the time index of the data read, also used in warning messages.
#' @param fill	is the value to use for missing data (outside of the mask).
#' @param keep.vxIX	is TRUE to return also the indices at the non-missing acoustic data.
#'
#' @importFrom TSD mergeListKeepDimensions
#'
#' @export
#' @rdname read.event_unzip_vbsc
#'
read.event_unzip_vbsc <- function(x, pad=TRUE, split=TRUE, filename=NULL, t=NULL, fill=NA, keep.vxIX=TRUE){
	# If no vbsc data har positive length, fill with 'fill' or do nothing if fill has length 0:
	if(length(x$vbsc) == 0){
		warn <- paste0("Volume backscattering coefficient 'vbsc' not present", 
				if(length(filename)) paste0(" in file ", filename), 
				if(length(t)) paste0(" for time ", paste(t, collapse=", "))
				)
		if(length(fill)){
			warn <- paste0(warn, "(set to an array of '", fill, "')")
			numt <- length(x$utim)
			x$vbsc <- vector("list", numt)
			for(i in seq_len(numt)){
				x$vbsc[[i]] = array(fill, dim=c(max(x$lenb[,i]), x$numb[i]))
			}
		}
		warning(warn)
	}
	
	# Filling in for compressed acoustic data:
	if(length(x$vbsc)>0 && length(x$vxIX)>0){
		# In the unlikely event that all time steps have equal number of positive voxels in the compressed mode, split into a list:
		if(!is.list(x$vbsc)){
			thislen <- nrow(x$vbsc)
			thisnumt = ncol(x$vbsc)
			splitvec <-  rep(seq_len(thisnumt), each=thislen)
			x$vbsc = split(x$vbsc, splitvec)
			x$vxIX = split(x$vxIX, splitvec)
			#x$vbsc = as.data.frame(x$vbsc)
			#x$vxIX = as.data.frame(x$vxIX)
			}
		# Create a list of vbsc, where each ping is represented in each element of the list, and fill inn the non-empty values indicated by 'vxIX':
		#for(i in seq_along(x$vbsc)){
		numt <- max(length(x$vbsc), length(x$utim))
		for(i in seq_len(numt)){
			#tempvbsc = NAs(max(x$lenb[,i]), x$numb[i])
			tempvbsc = array(fill, dim=c(max(x$lenb[,i]), x$numb[i]))
			# Use c() since 'lenb' and 'numb' are required to be arrays with time along last dimension, obtained using drop.out=FALSE in read.TSD()
			tempvbsc[c(x$vxIX[[i]])] = c(x$vbsc[[i]])
			x$vbsc[[i]] = tempvbsc
		}
		# Remove vxIX if present:
		if(!keep.vxIX){
			x$vxIX <- NULL
		}
	}
	
	# Rearranging sv-values in a 3 dimensional array [lenb, numb, numt]:
	if(is.list(x$vbsc)){
		x$vbsc = lapply(x$vbsc, function(x) if(length(dim(x)) == 2) array(x, dim=c(dim(x),1)) else x)
		}
	else{
		lenb <- x$lenb[1]
		numb <- x$numb[1]
		dim(x$vbsc) <- c(lenb, numb, length(x$vbsc) / (lenb*numb))
	}
	# Check whether the dimensions of each time step are identical:
	if(is.list(x$vbsc)){
		x$vbsc = mergeListKeepDimensions(x$vbsc, pad=pad, split=split, add1=length(dim(x$vbsc[[1]])) == 2)
		}
	return(x)
}
