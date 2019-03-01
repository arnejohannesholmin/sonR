#*********************************************
#*********************************************
#' (Internal) Returns a summary of segmentation data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ang2rot combine.TSD zeros
#' @importFrom stats quantile
#'
#' @export
#' @rdname summarySeg
#' 
summarySeg = function(event, t, TS=-55, w=0.369, segpar=1:2, ang=NULL, maxRange=500, n90=2){
	outlist = vector("list", 2)
	out = zeros(3, 12)
	for(i in seq_along(segpar)){
		s = read.event(event, t=t, var="seg", segpar=segpar[[i]])
		outlist[[i]] = combine.TSD(s)
		
		get.outliers = function(x, n90){
			u90 = quantile(x, 0.9, na.rm=TRUE)*n90
			cat("Upper outlier threshold: ", u90, "\n")
			(x > u90) | (x < 0)
			}
		
		# Add angle, tons, tilt
		outlist[[i]]$tons = s$Xtbs / 10^(TS/10) * w * 1e-3
		outlist[[i]]$tonsE = s$Xebt / 10^(TS/10) * w * 1e-3
		outlist[[i]]$anio = ang2rot(s$anio)*180/pi
		outlist[[i]]$anis = ang2rot(s$anis)*180/pi
		outlist[[i]]$tilt = 180 - read.event(event, var="dire", t=t)$dire[1,]*180/pi
		
		
		# Remove extreme and negative points:
		outliers = get.outliers(outlist[[i]]$tons, n90) | is.na(outlist[[i]]$tons)
		
		presentvar = c("tons", "tonsE", "Xtbs", "XtTS", "Xasv", "XaSv", "Xtvl", "Xtha")
		for(j in seq_along(presentvar)){
			outlist[[i]][[presentvar[j]]][outliers] = NA
			}
		
		if(i==1){
			if(length(ang)==1){
				ang = c(ang, ang + 180 * floor(max(outlist[[i]]$anio-ang[1])/180))
				}
			ang_subset = seq_along(outlist[[i]]$tons)
			if(length(ang)>1){
				ang_subset = which(outlist[[i]]$anio>=ang[1] & outlist[[i]]$anio<=ang[2])
				}
			range_subset = seq_along(outlist[[i]]$tons)
			if(length(maxRange)>0){
				range_subset = which(outlist[[i]]$Xhra < maxRange)
				}
			subset = intersect(ang_subset, range_subset)
			}
		
		
		out[1,0+i] = mean(outlist[[i]]$tons[subset], na.rm=TRUE)
		out[2,0+i] = quantile(outlist[[i]]$tons[subset], 0.5, na.rm=TRUE)
		out[3,0+i] = quantile(outlist[[i]]$tons[subset], 0.9, na.rm=TRUE)
		
		out[1,2+i] = 10*log10(mean(outlist[[i]]$Xtbs[subset], na.rm=TRUE))
		out[2,2+i] = quantile(outlist[[i]]$XtTS[subset], 0.5, na.rm=TRUE)
		out[3,2+i] = quantile(outlist[[i]]$XtTS[subset], 0.9, na.rm=TRUE)

		out[1,4+i] = 10*log10(mean(outlist[[i]]$Xasv[subset], na.rm=TRUE))
		out[2,4+i] = quantile(outlist[[i]]$XaSv[subset], 0.5, na.rm=TRUE)
		out[3,4+i] = quantile(outlist[[i]]$XaSv[subset], 0.9, na.rm=TRUE)

		out[1,6+i] = mean(outlist[[i]]$Xtvl[subset], na.rm=TRUE)
		out[2,6+i] = quantile(outlist[[i]]$Xtvl[subset], 0.5, na.rm=TRUE)
		out[3,6+i] = quantile(outlist[[i]]$Xtvl[subset], 0.9, na.rm=TRUE)

		out[1,8+i] = mean(outlist[[i]]$Xtha[subset], na.rm=TRUE)
		out[2,8+i] = quantile(outlist[[i]]$Xtha[subset], 0.5, na.rm=TRUE)
		out[3,8+i] = quantile(outlist[[i]]$Xtha[subset], 0.9, na.rm=TRUE)

		out[1,10+i] = mean(outlist[[i]]$tonsE[subset], na.rm=TRUE)
		out[2,10+i] = quantile(outlist[[i]]$tonsE[subset], 0.5, na.rm=TRUE)
		out[3,10+i] = quantile(outlist[[i]]$tonsE[subset], 0.9, na.rm=TRUE)
		}
	out = rbind(out, out[3,] / out[2,])
	colnames(out) = c("W1", "W2", "TS1", "TS2", "Sv1", "Sv2", "V1", "V2", "A1", "A2", "WE1", "WE2")
	rownames(out) = c("Mean", "50%-ile", "90%-ile", "Ratio_90_50")
	out = round(out, digits=2)
	outlist$tab = out
	outlist$ang = ang
	outlist
	}
