#*********************************************
#*********************************************
#' (Internal) Get midpoints of voxels.
#'
#' @param x		A list of data as the output from read.event(), containing the variables 'lenb', 'dire', 'psxv', 'psyv', 'pszv', and those required by \code{\link{soundbeam.TSD}}.
#' @param pad	Logical: If TRUE, pad with NAs for non-equal lengths of beams.
#' @param split	See \code{\link{mergeListKeepDimensions}}.
#' @param ...	Passed on to \code{\link{soundbeam.TSD}}.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD dim_all mergeListKeepDimensions numt.TSD
#'
#' @export
#' @rdname psx.TSD
#' 
psx.TSD <- function(x, pad=TRUE, split=TRUE, ...){
	
	if(isTRUE(pad)){
		x$lenb <- array(max(x$lenb), dim=dim_all(x$lenb))
	}
	# If vessel positions are missing, set these to 0:
	numt <- numt.TSD(x)
	if(length(x$psxv)==0){
		x$psxv <- double(numt)
		x$psyv <- double(numt)
		x$pszv <- double(numt)
	}
	
	out <- rep(list(list()),3)
	names(out) <- c("psxx", "psyx", "pszx")
	# Midpoints of voxels:
	# For echosounders speed up by only applying soundbeam.TSD() once, and reapeating for all time steps, added the vessel positions:
	tol <- 0.05
	if(pad && all(abs(x$dire-pi)<tol)){
		# Get voxel positions for the first ping:
		temp <- soundbeam.TSD(data=x, t=1, plot=FALSE, rpos="midp", ...)
		dimout <- c(dim_all(temp$x), numt)
		lengthOne <- prod(dim_all(temp$x))
		# Repeat for the other pings, adding vessel positions (restrict to valid pings, if the length if the positions does not match the mtim, a problem that may appear if position data is only available via adds in plotting functions, and we are moving beyond the time range of the data):
		valid <- seq_len(numt)
		out$psxx <- c(temp$x) + rep(x$psxv[valid], each=lengthOne) - x$psxv[1]
		out$psyx <- c(temp$y) + rep(x$psyv[valid], each=lengthOne) - x$psyv[1]
		out$pszx <- c(temp$z) + rep(x$pszv[valid], each=lengthOne) - x$pszv[1]
		dim(out$psxx) <- dimout
		dim(out$psyx) <- dimout
		dim(out$pszx) <- dimout
		return(out)
	}
	else{
		for(i in seq_len(numt)){
			temp <- soundbeam.TSD(data=x, t=i, plot=FALSE, rpos="midp", ...)
			out$psxx[[i]] <- temp$x
			out$psyx[[i]] <- temp$y
			out$pszx[[i]] <- temp$z
		}
		return(lapply(out, mergeListKeepDimensions, pad=pad, split=split, add1=length(dim(out$psxx[[1]]))==2))
	}
}
