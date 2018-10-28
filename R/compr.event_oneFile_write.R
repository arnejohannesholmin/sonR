#*********************************************
#*********************************************
#' (Internal) Read and write one Simrad raw file.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD write.TSD numt.TSD
#'
#' @export
#' @rdname EKRaw2TSD_oneFile_write
#' 
compr.event_oneFile_write <- function(i, indt, filelist, pingsfiles, vesselfiles, beamsfiles, t="all", drop=TRUE, compress=FALSE, TIME=NULL, write=TRUE, 
	# Inputs used in compr.TSD:
	tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("mean","median"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, z0=0, pingsnames=NULL, vesselnames=NULL, beamsnames=NULL, dumpfiles=NULL, maxlenb=NULL, ...){

	# Function used for writing the first, last and the rest of the time steps in three separate files:
	write.TSD123 <- function(data, names, files, i, numt, ...){
		# Write the first ping:
		if(length(dumpfiles)){
			write(paste("TSD1_", i), dumpfiles[i], append=TRUE)
		}
		this = subset_TSD(data[names], ind=list(1), pad="start", drop=FALSE)
		TSD::write.TSD(this, files$comprFileNames1[i], numt=1, ...)
		# Write the middle pings:
		if(numt>1){
			if(length(dumpfiles)){
				write(paste("TSD2_", i), dumpfiles[i], append=TRUE)
			}
			thisindt <- seq(2, max(2,numt-1))
			this = subset_TSD(data[names], ind=list(thisindt), pad="start", drop=FALSE)
			TSD::write.TSD(this, files$comprFileNames2[i], numt=length(thisindt), ...)
		}
		# Write the last ping:
		if(numt>2){
			if(length(dumpfiles)){
				write(paste("TSD3_", i), dumpfiles[i], append=TRUE)
			}
			this = subset_TSD(data[names], ind=list(numt), pad="start", drop=FALSE)
			TSD::write.TSD(this, files$comprFileNames3[i], numt=1, ...)
		}
	}
		
	# Read the data from the raw file:
	#cat("Reading pings file number ", i, " of ", length(pingsfiles$filenames), " (", basename(pingsfiles$filenames[i]), ") ...\n", sep="")
	if(length(dumpfiles)){
		write(paste("File", i), dumpfiles[i], append=TRUE)
	}
		
	#data = read.event(filelist, t=indt[[i]], var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE, onestep=FALSE)
	data = read.event(filelist, t=indt[[i]], var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE, onestep=2)
	
	# Discard pings with too long beams:
	maxLenbPerPing = if(length(dim(data$lenb))==2) apply(data$lenb, 2, max) else max(data$lenb)
	if(length(maxlenb)>0 && any(maxLenbPerPing>maxlenb)){
		validPings = maxLenbPerPing <= maxlenb
		data = read.event(filelist, t=indt[[i]][validPings], var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE)
		write(paste("Discarded pings:", paste(which(!validPings), collapse=", ")), dumpfiles[i], append=TRUE)
	}
	#data = read.event(filelist, t=head(indt[[i]], 3), var=c(pingsnames, beamsnames, vesselnames), TIME=TIME, other=FALSE, drop.out=FALSE)
	
	# Abort empry data:
	numt = numt.TSD(data)
	if(numt==0){
		warning(paste("No valid pings in file", i))
		return()
	}
	# Compress the data:
	if(compress){
		if(length(dumpfiles)){
			write(paste("Compr", i), dumpfiles[i], append=TRUE)
		}
		data <- compr.TSD(data, tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, z0=z0, ...)
		# Update numt for the compressed data:
		numt = numt.TSD(data)
	}
	
	if(length(dumpfiles)){
		write(paste("Write", i), dumpfiles[i], append=TRUE)
	}
	if(write){
		# Write the first time step, the last time step and the center time steps, so that the last and first time steps can be merged between data coming from consecutive raw files:
		if(compress){
			write.TSD123(data, pingsnames, pingsfiles, i, numt, header=list(dtyp=list(vbsc="floa")))
			write.TSD123(data, beamsnames, beamsfiles, i, numt)
			write.TSD123(data, vesselnames, vesselfiles, i, numt)
		}
		else{
			TSD::write.TSD(data[pingsnames], pingsfiles$filenames[i], numt=numt, header=list(dtyp=list(vbsc="floa")))
			TSD::write.TSD(data[beamsnames], beamsfiles$filenames[i], numt=numt)
			TSD::write.TSD(data[vesselnames], vesselfiles$filenames[i], numt=numt)
		}
	}
	if(length(dumpfiles)){
		write(paste("End", i), dumpfiles[i], append=TRUE)
	}
	#invisible(data)
	rm(data)
	gc()
	return(filelist[i])
}
