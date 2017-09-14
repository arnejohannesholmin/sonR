#*********************************************
#*********************************************
#' (Internal) Read and write one Simrad raw file.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD labl.TSD write.TSD
#'
#' @export
#' @rdname EKRaw2TSD_oneFile_write
#' 
EKRaw2TSD_oneFile_write <- function(i, filelist, pingsfiles, vesselfiles, rawvesselfiles, beamsfiles, pingsfiles1=NULL, pingsfiles2=NULL, pingsfiles3=NULL, prenumt=10, t="all", endian="little", timeOffset=0, minTimeDiff=Inf, msg=TRUE, na.rm=TRUE, correctTime=TRUE, bmmd=NULL, TVG.exp=2, dira_offset=0, compress=FALSE, write=TRUE, cali=TRUE, toTS=FALSE, psze=-7, cleanNMEA=1, 
	# Inputs used in compr.TSD:
	tres=NULL, xres=NULL, zres=NULL, rres=NULL, bres=NULL, funvbsc=c("median","mean"), funt=c("median","mean"), adds=list(), split=TRUE, skipAngles=TRUE, origin=1, apply.range.offset=FALSE, thr1m=FALSE, ...){

	# Declare the variable names:
	beamsnames = TSD::labl.TSD("EKRaw2TSD_b")
	vesselnames = TSD::labl.TSD("EKRaw2TSD_v")
	rawvesselnames = TSD::labl.TSD("EKRaw2TSD_r")
	pingsnames = TSD::labl.TSD("EKRaw2TSD_p")

	# Read the data from the raw file:
	data = EKRaw2TSD_oneFile(i=i, filelist=filelist,  prenumt=prenumt, t=t, endian=endian, timeOffset=timeOffset, minTimeDiff=minTimeDiff, msg=msg, na.rm=na.rm, correctTime=correctTime, bmmd=bmmd, TVG.exp=TVG.exp, dira_offset=dira_offset, cali=cali, toTS=toTS, psze=psze, skipAngles=skipAngles, cleanNMEA=cleanNMEA, apply.range.offset=apply.range.offset, thr1m=thr1m)
	if(length(data)==0){
		rm(data)
		gc()
		return(i)
		}
	numt = length(data$mtim)
		
	# Compress the data:
	if(compress){
		data <- compr.TSD(data, tres=tres, xres=xres, zres=zres, rres=rres, bres=bres, funvbsc=funvbsc, funt=funt, adds=adds, split=split, skipAngles=skipAngles, origin=origin, ...)
		}
	
	# Write the first time step, the last time step and the center time steps, so that the last and first time steps can be merged between data coming from consecutive raw files:
	if("p" %in% write){ 
		if(compress){
			# Write the first ping:
			TSD::write.TSD(data[pingsnames], pingsfiles1[i], t=1, numt=numt, header=list(dtyp=list(vbsc="floa")))
			# Write the middle pings:
			TSD::write.TSD(data[pingsnames], pingsfiles2[i], t=seq(2,numt-1), numt=numt, header=list(dtyp=list(vbsc="floa")))
			# Write the last ping:
			TSD::write.TSD(data[pingsnames], pingsfiles3[i], t=numt, numt=numt, header=list(dtyp=list(vbsc="floa")))
			}
		else{
			TSD::write.TSD(data[pingsnames], pingsfiles[i], numt=numt, header=list(dtyp=list(vbsc="floa")))
			}
		}
	
	# Write the remaining data regardless of compression of the acoustic data:
	if("b" %in% write){ 
		TSD::write.TSD(data[beamsnames], beamsfiles[i], numt=numt)
		}
	if("v" %in% write){ 
		TSD::write.TSD(data[vesselnames], vesselfiles[i], numt=numt, header=list(dtyp=list(mtim="doub", lonv="doub", latv="doub", sadv="doub")))
		}
	if("rv" %in% write){ 
		TSD::write.TSD(data[rawvesselnames], rawvesselfiles[i], numt=1, header=list(dtyp=list(imtm="doub", ilnv="doub", iltv="doub", isdv="doub")))
		}

	rm(data)
	gc()
	#invisible(data)
	i
	}
