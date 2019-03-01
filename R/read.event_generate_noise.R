#*********************************************
#*********************************************
#' Estimates the noise in MS70 data.
#'
#' @param data  is a list containing the data from which the noise should be estimated. Sould contain the following variables: "vbsc", "time", "hini", "hins", "beams", "ctd".
#' @param t  is either the indexes of the pings to be treated, as listed from 1 to the number of pings in the event, or the time point given as a string "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF".
#' @param noise  is a vector of strings representing the noise components to include in the estimate.
#' @param nsind  is a vector of indexes along the beams, as input to ind.expand(), used to select the subset over which the estimation of the phase of the periodic noise is done. If given as a single numeric, the outermost 'nsind' voxels are used in each beam.
#' @param cruise  Used by \code{\link{get.pdns_phase.TSD}} to get the phase of the periodic noise.
#' @param hins_add  is the number of voxels that should be discarded on both sides of high intensity noise voxels voxels along beams, used for accounting for possible high values that are related to the high intensity noise but not classified as such voxels.
#' @param phase,esnm,dir.data  is FALSE if any of 'pn3M' (phase for each time step) or 'pns3' (phase equal for all time steps) given in 'bgns' or read from the noise file located by the funciton noise.path.event() should be used, as oposed to estimating the phase from the data for each time step. This is only recommended for simulated data where the phase is constant over all time steps, and saves some CPU time.
#' @param pdns_scale  is used in get.pdns_phase.event() to scale the noise in order to allow the optimization to work.
#' @param TVG  is FALSE if TVG compensation is to be removed from the data.
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw apply.TVG
#' @importFrom TSD arr.ind2ind zeros
#'
#' @export
#' @rdname read.event_generate_noise
#'
read.event_generate_noise<-function(data, t=1, noise=c("bgns","pdns","nrns","hins"), nsind=0.75, cruise=2009116, esnm="MS70", dir.data=NULL, hins_add=10, phase=TRUE, pdns_scale=1e-14, TVG=TRUE, TVG.exp=2){
	
	############### LOG: ###############
	# Start: 2012-11-23 - Clean version, adopted from echoIBM.vbsc2p.event().
	
	##### Preparation #####
	# Store the original dimension of the acoustic data:
	if(length(dim(data$vbsc))==2){
		dim(data$vbsc) = c(dim(data$vbsc),1)
		}
	olddim=dim(data$vbsc)
	
	# Define the output total noise:
	data$tlns = zeros(olddim)
	namesdata = names(data)
	
	bgnsPresent = any(c("bgns", "bgn0") %in% namesdata)
	pdnsPresent = any(c("pns1", "pn3M", "pns3") %in% namesdata)
	nrnsPresent = any(c("nr0a", "nr0p") %in% namesdata)
	hinsPresent = length(data$hini)>0
	
	# Expand the background noise to an array of the same size as the data:
	if("bgns" %in% noise && bgnsPresent){
		thisbgns = data[c("bgns","bgn0")]
		thisbgns = thisbgns[[which(sapply(thisbgns, function(x) length(x)>0))[1]]]
		data$tlns = data$tlns + c(matrix(thisbgns, nrow=max(data$lenb),ncol=data$numb,byrow=TRUE))
		}
	
	# If the required information is present in 'data', estimate the phase of the periodic noise at each time step, from a portion of the sonar volume dominated by noise, or alternatively, not occupied by schools:
	if("pdns" %in% noise && pdnsPresent){
		periodic = which(data$badb==1)
		
		if(phase){
			if(length(data$acfq)>0 && !any(is.null(data$pns1),is.null(data$pns2),is.null(data$harm))){
				data$pn3M = get.pdns_phase.TSD(data=data, nsind=nsind, cruise=cruise, esnm=esnm, dir.data=dir.data, pdns_scale=pdns_scale, TVG=TVG)$pn3M
				}
			else{
				stop("Estimating the phase of the periodic noise requires 'acfq', 'pns1', 'pns2' and 'harm' to be present in 'bgns' or read from the noise files located in the central noise directory of the cruise given by 'cruise'")
				}
			}
		else{
			if(ncol(data$pn3M)<max(t)){
				warning(paste("The number of columns in 'pn3M' (",ncol(data$pn3M),") should be at least equal to the number of pings, Trying 'pns3'..."))
				if(length(data$pns3)>0){
					data$pn3M = matrix(data$pns3,nrow=length(data$freq),ncol=max(t))
					}
				else{
					stop("The phase of the periodic noise was not found in the data. Try phase=TRUE to estimate the phase at each ping")
					}
				}
			}
		}
	
	# If the required information is present in 'data', add near range noise:
	if("nrns" %in% noise && nrnsPresent){
		if(length(data$nr0a)>0){
			data$tlns = data$tlns + c(data$nr0a[nrow(data$tlns), ])
			}
		else if(length(data$nr0p)>0){
			data$tlns = data$tlns + c(data$nr0p[nrow(data$tlns), ])
			}
		}
	
	
	##### Execution and output #####
	# Move through the time steps:
	if("pdns" %in% noise && pdnsPresent){
		for(i in seq_along(t)){
			if("pdns"%in%noise){
				# Estimate the expected periodic noise:
				data$pdns = echoIBM.pdns(data,indt=i)$pdns
				# Add the periodic noise:
				data$tlns[,periodic,i] = data$tlns[,periodic,i] + data$pdns[,periodic]
				}
			}
		}
		
	# Add TVG to the total noise thus far, before adding the high intensity noise, which is taken from the TVG-amplified data:
	data$tlns = apply.TVG(data$tlns,data,TVG.exp=TVG.exp)
		
	
	# Add the high intensity noise:
	if("hins"%in%noise && hinsPresent){
		for(i in seq_along(t)){
			# Select only the high intensity noise of the current time step, and skip if the high intensity noise is not present:
			thesehini = data$hini[data$hini[,3]==t[i],1:2,drop=FALSE]
			if(length(thesehini)>0){
				# Expand the high intensity noise along the beams:
				if(hins_add>0){
					
					newhini = NULL
					uniquebeams = unique(thesehini[,2])
					for(b in uniquebeams){
						thesej = sort(thesehini[thesehini[,2]==b,1])
						diffj = diff(thesej)
						jump = which(diffj>1)
						if(length(jump)==0){
							lower = 1
							upper = length(thesej)
							}
						else{
							lower = c(1,jump+1)
							upper = c(jump,length(thesej))
							}
						# Expand the high intensity noise voxels along beams:
						expanded = NULL
						for(l in seq_along(lower)){
							expanded = c(expanded, (thesej[lower[l]]-hins_add):(thesej[upper[l]]+hins_add))
							}
						# Uniquify the expanded voxels:
						expanded = intersect(expanded,seq_len(max(data$lenb)))
						expanded = cbind(unique(expanded),b)
						newhini = rbind(newhini,expanded)
						}
					thesehini = cbind(newhini,i)
					}
				
				thesehini = arr.ind2ind(thesehini,c(max(data$lenb),data$numb,length(t)))
				data$tlns[thesehini] = data$vbsc[thesehini]
				}
			}
		}
		
	# Return the output:
	data$tlns
}
