#*********************************************
#*********************************************
#' Estimates and writes the noise in av underwater acoustic system like Simrad EK60, ME70, MS70 or SX90, based on Holmin et al. 2013: "Estimation and simulation of multibeam sonar noise".
#'
#' @param con  is the path to the file to which the noise estimates schould be written. If con==NULL, no data is written to file (only returned from the function). If con==TRUE, the data is written to the default file located in the directory "Noise/Main" in the cruise directory, and named with the cruise name and the events used to estimate the various noise types.
#' @param cruise  is the identifier of the cruise for the passive noise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. Currently implemented are "EK60", "ME70", "MS70" and "SH80"/"SX80"/"SH90"/"SX90" (may be given in lover case).
#' @param treat  is a vector of strings naming the noise types to estimate, where "bgns" is the background noise, "pdns" is the periodic noise, "nrns" is the near range noise, and "hins" is the high intensity noise.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param event_p  is the identifier of the event for the noise in passive mode. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
#' @param t_p  is the identifier of the time points for the noise in passive mode. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
#' @param noise  is either a list of the elements 'harm' and 'badb', representing the harmonics of the periodic beams (where the non-periodic beams are set to 1) and the the vector identifying the periodic beams, composed of 0's for non-periodic beams and 1's for periodic beams, respectively; or a character string representing the path to a directory or file from which this information should be read. If not given, the periodic nature and the harmonics of the beams are estimated from the passive data.
#' @param fftthreshold  is the threshold value for the height of the highest peak in spectral density relative to the zero frequency, identifying the preiodic beams as those with mean across pings exceeding the 'fftthreshold'.
#' @param nsind  is a specification of the voxels used in the estimation, given either as a list as input to pplot3d.event(), a vector to be intersected with the sequence of voxels along beams, a single integer > 1 representing the "farthes nsind voxels along each beam", a single numeric <=1 representing the proportion of the farthest voxels along the beams compared to the length of the beams, of sonething else (like NULL) implying taking all voxels into account:.
#' @param acfq  is the frequency of the alternating current generating the periodic noise.
#' @param prob  is the lower fraction of tha data in each beam to be used in the estimation of background noise. Lower value gives higher variance in the estimate and higher value gives stronger influence from periodic noise.
#' @param type  is the type of quantile method used in meanexp.quantile() in get.bgns.event(). Defaulted to 6, which was found reasonable for exponential data.
#' @param tries  specifies the number of tries in the iteration used when estimating the background noise and the periodic noise (vector of two values where the first limits the iteration between estimation of the background noise and the periodic noise and the second limits the oprimization of the periodic noise).
#' @param thr  is the threshold for the absoulte difference between the previous and the current estimate, divided by the current estimate, at which the above mentioned iteration stops.
#' @param trim  is used in mean() for estimation of the kurtosis parameter of the periodic noise 'pns2'.
#' @param search  is the number of tries in the iteration above which the optimal estimate is searched for in the intermediate estimates, by the value with the lowest sum of the rank of the funciton value of the iteration (bgns) and of the optimization (pdns).
#' @param kurtosis  is a vector of the lower and upper limit for the kurtosis estimate, outside of which estimates (of all parameters) are discarded.
#' @param reffan  is the reference fan used in the calculation of the phase angles.
#' @param scale  is used for scaling the data to ensure best performance of the optimization, and corersponds to the mean of the data (default if set to NULL).
#' @param lowerpar , 'startpar' and 'upperpar' are the lower, start and upper parameter values representing the magnitude, kurtosis and phase, respectively, where the magnitude is given relative to 'scale'.
#' @param kern  is either a vector of wieghts summing to 1, specifying the smoothing of the spectral frequency used to identify the harmonics, or a single numeric = 1, specitying a Gaussian kernel with sd=1.
#' @param high  is used in resampleHigh() in get.bgns.event() to reduce the impact from unusually high values in the fft.
#' @param event_a  is the identifier of the event for the noise in active mode. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
#' @param t_a  is the identifier of the time points for the noise in active mode. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
#' @param filter  is the kernel used when smoothing the near range noise of each beam.
#' @param nsdn  is the number of standard deviation above the mean nrea range noise for each beam, defining the values that are set to 0. The beam is traced from the sonar and out, and at the first value below the threshold, the near range noise is set to zero for all values at and beyond this value.
#' @param farv  is the number of far voxels used for the calclulation of the mean and standard deviation used when assigning zeros to the smoothed near range noise. For farv=13, the 13 farthest non-missing values are used for each beam.
#' @param surface  is a vector of indexes for the beams affected by the surface noise (extrapolated in the near range noise estimation).
#' @param max.memory  puts a restiction on the mempry occupied by the function, prompting a call to the user to approve is exceeded.
#' @param turns  is a numeric giving the length of the runs, causing 'turns' pings to be read at each step in the function.
#' @param k  is a numeric giving the width of the medians across pings, used in runmed() inside get.hins.event_small().
#' @param q  is either a single numeric, or a vector of length 2, givin the start and end point in the linear function along beams defining the threshold times the median filtered values in each direction above which data are classified as high intensity noise.
#' @param beta_school  is a single numeric representing a typical high school value. Only values above 'beta_school' can be classified as high intensity noise. 'beta_school' assures that spikes due to for example fish in very silent regions of the data are not classified as high intensity noise.
#' @param ind  is a list of indexes, as given to subset_TSD(), used to select the subset over which the estimation of high intensity noise is done. Defaulted to exclude the first 100 voxels along each beam.
#' @param ...  variables used in write.TSD() (such as 'ow' for overwriting existing files).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD labl.TSD write.TSD
#' @importFrom utils tail
#' @importFrom stats dnorm
#'
#' @export
#' @rdname write.noise.cruise
#'
write.noise.cruise<-function(
	# Main variables:
	con=NULL, cruise=NULL, esnm=NULL, treat=c("bgns","pdns","nrns","hins"), dir.data=NULL, 
	# Used in get.bgns.event():
	event_p=NULL, t_p=NULL, noise=NULL, fftthreshold=0.145, nsind=0.75, acfq=100, prob=0.2, type=6, tries=c(10,10), thr=1e-3, trim=0.1, search=4, kurtosis=c(0.1,3), reffan=16, scale=1e-14, lowerpar=c(1e-3,0.001,-2*pi), startpar=c(1,2,3), upperpar=c(1e3,6,4*pi), kern=NULL, high=Inf, 
	# Used in get.nrns.event():
	event_a=NULL, t_a=NULL, filter=dnorm(-20:20,sd=4), nsdn=1000, farv=1000, surface=451:500, 
	# Used in get.hins.event():
	turns=114, k=5, q=c(1e3,1e2), beta_school=2e-4, ind=list(-(1:100),NULL), max.memory=2e9, 
	# Passed on to write.TSD():
	...){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-03-13 - Clean version.
	# Update: 2010-06-03 - Changed to be able to estimate noise from several events.
	# Update: 2012-05-23 - Updated to accept information about the periodic noise.
	# Update: 2012-05-23 - Updated to use the newest version of get.bgns.event, which fistly identifies the beams affected by periodic noise, secondly estimates the background noise by simple mean in the healthy beams and conditional mean in the beams affected by the periodic noise, and at last estimates the parameters of the periodic noise modeled by a1 exp(a2 sin(n * 2pi * f * j - a3)).
	# Update: 2012-06-15 - Changed method to specifying passive and active noise events, and calclulate the background noise from the passive event, and near range noise from both types of events. Also more detailed output given.
	# Update: 2012-06-29 - Changed to use the info variable returned from get.bgns().
	# Update: 2012-10-26 - Removed defaults for 't', 'event' and 'cruise'.
	# Update: 2012-10-29 - Added the option 'lowerkurtosis', which discards the estimations for which the kurtosis estimate is lower than 0.5, to avoid erroneous estimates from affecting the mean estimates.
	# Update: 2012-10-31 - Added the options 'scale', 'lowerpar', 'startpar' and 'upperpar'.
	# Update: 2012-11-19 - Changed the name of 'lowerkurtosis' to 'kurtosis', and added the option of both lower and upper limit for the kurtosis in this parameter.
	# Update: 2012-11-24 - Added the simple mean of each beam, regardless of periodic noise, to the output.
	# Update: 2012-12-31 - Removed writing info, since the info is added when the data are read.
	# Last: 2013-07-19 - Removed 'near.range' and 'ind' from get.nrns.event() and renamed 'ind1' to 'ind' for use in get.bgns.event().
		
	
	##################################################
	##################################################
	########## Preparation ##########
	### if(length(con)==0){
	### 	stop("'con' must be given")
	### 	}
	# If "pdns" is not present in treat, set harm to 'ones' and 'badb' to zeros:
	if(!"pdns" %in% treat){
		harm = 1
		badb = 0
		}
	# Define the names of the variables to write/return:
	bgnsnames = labl.TSD("bgns")
	nrnsnames = labl.TSD("nrns")
	nrnpnames = labl.TSD("nrnp")
	nrnanames = labl.TSD("nrna")
	hinsnames = labl.TSD("hins")
	# Background noise is needed for estimation of the near range noise:
	if("nrns" %in% treat && !"bgns" %in% treat){
		treat = c("bgns", treat)
		}
		
	# Get the path to the cruise:
	cruise = event.path(cruise=cruise, dir.data=dir.data, ...)

	
	########## Execution ##########
	# Estimate the background noise using the echoIBM utility function get.bgns.event(), for all pings in the time range 't_bgns':
	if("bgns" %in% treat){
		bgns = get.bgns.event(con=NULL, event=event_p, cruise=cruise$cruise, t=t_p, esnm=esnm, noise=noise, fftthreshold=fftthreshold, nsind=nsind, acfq=acfq, prob=prob, type=type, tries=tries, thr=thr, trim=trim, search=search, kurtosis=kurtosis, reffan=reffan, scale=scale, lowerpar=lowerpar, startpar=startpar, upperpar=upperpar, kern=kern, high=high)
		# Select only the bgns-variables to write/return:
		if(length(bgns)>0){
			bgns = bgns[bgnsnames]
			}
		if(any(is.na(bgns$bgns))){
			warning("The kurtosis of the periodic noise was not correctly estimated (below fftthreshold = ", fftthreshold, ") for the following beams:", paste0(which(is.na(bgns$bgns)), collapse=", "))
			}
		}
	else{
		bgns = list()
		}
	
	# Estimate the near range noise for the passive data, without surface compensation:
	if("nrns" %in% treat){
		nrnp = get.nrns.event(con=NULL, event=event_p, cruise=cruise$cruise, t=t_p, bgns=bgns, filter=filter, nsdn=nsdn, farv=farv, surface=NULL, esnm=esnm)
		# If the same pings are selected for active and passive mode, copy the results from above:
		if(identical(event_p, event_a) && identical(t_p, t_a)){
			nrna = nrnp
			}
		else{
			# Estimate the near range noise for the active data, with surface compensation:
			nrna = get.nrns.event(con=NULL, event=event_a, cruise=cruise$cruise, t=t_a, bgns=bgns, filter=filter, nsdn=nsdn, farv=farv, surface=surface, esnm=esnm)
			}
		
		# Select only the nrnp-variables to write/return:
		if(length(nrnp)>0){
			nrnp = nrnp[nrnsnames]
			if(length(nrnp)>0){
				# Change names to include the mode of the data:
				names(nrnp) = nrnpnames
				}
			}
		
		# Select only the nrna-variables to write/return:
		if(length(nrna)>0){
			nrna = nrna[nrnsnames]
			if(length(nrna)>0){
				# Change names to include the mode of the data:
				names(nrna) = nrnanames
				}
			}
		}
	else{
		nrnp = list()
		nrna = list()
		}
	# Estimate the high intensity noise:
	if("hins" %in% treat){
		# Estimate the high intensity noise for all events:
		events = list.files(file.path(cruise$cruise, "Events"), full.names=TRUE)
		hinsfilenames = substring(events, sapply(gregexpr("/", events), tail, 1)+1)
		events = file.path(events, esnm, "tsd")
		hinsfilenames = file.path(events, paste0(hinsfilenames, "_hins.beams"))
		# Run through the events and detect the high intensity noise:
		#HINS = vector("list",length(hinsfilenames))
		#HINS = rep(list(HINS), 6)
		#names(HINS) = c("hins", "hini_j", "hini_i", "hini_indt", "hini_eventnr", "hini_mtim")
		for(i in seq_along(hinsfilenames)){
			hins = get.hins.event(con=hinsfilenames[i], event=events[i], t="all", turns=turns, k=k, q=q, beta_school=beta_school, ind=ind, max.memory=max.memory)
			#nhins = length(thishins$hins)
			#if(nhins>0){
			#	# Collapse the data into the matrix HINS, holding in the first column, the estimated high intensity noise, in the second through the fifth column, the along beam indices (j), the beam indices (i), the time step indices (indt) and the event indices, and in the sixth column the MATLAB serial date number (mtim):
			#	HINS = list(HINS=cbind(thishins$hins, thishins$hini[,1], thishins$hini[,2], thishins$hini[,3], rep(i,nhins), thishins$mtim[thishins$hini[,3]]))
			#	write(HINS, )
			#}
			
			#if(nhins>0){
			#	HINS[["hins"]][[i]] = thishins$hins
			#	HINS[["hini_j"]][[i]] = thishins$hini[,1]
			#	HINS[["hini_i"]][[i]] = thishins$hini[,2]
			#	HINS[["hini_indt"]][[i]] = thishins$hini[,3]
			#	HINS[["hini_eventnr"]][[i]] = thishins$mtim[thishins$hini[,3]]
			#	HINS[["hini_mtim"]][[i]] = rep(i,nhins)
			#	#gc()
			#	}
			}
		# Collapse the data into the matrix HINS, holding in the first column, the estimated high intensity noise, in the second through the fifth column, the along beam indices (j), the beam indices (i), the time step indices (indt) and the event indices, and in the sixth column the MATLAB serial date number (mtim):
		#hins = list(HINS=cbind(unlist(HINS$hins),unlist(HINS$hini_j),unlist(HINS$hini_i),unlist(HINS$hini_indt),unlist(HINS$hini_eventnr),unlist(HINS$hini_mtim)))
		}
	else{
		hins = list()
		}
	
	
	########## Output ##########
	if(isTRUE(con)){
		event_p = as.numeric(event.path(event=event_p, cruise=cruise$cruise, esnm=esnm)$eventnr)
		event_a = as.numeric(event.path(event=event_a, cruise=cruise$cruise, esnm=esnm)$eventnr)
		con = file.path(cruise$cruise, "Noise", "Main")
		if(!file.exists(con)){
			suppressWarnings(dir.create(con, recursive=TRUE))
			}
		filename = ""
		if(length(event_p)>0 && "bgns" %in% treat){
			filename = paste0(filename, "_bgnsEvent_", event_p,  "_pdnsEvent_", event_p)
			}
		if(length(event_a)>0 && "nrns" %in% treat){
			filename = paste0(filename, "_nrnsEvent_", event_a)
			}
		if("hins" %in% treat){
			filename = paste0(filename, "_hinsAll")
			}
		con = file.path(con, paste0(substr(cruise$cruisename, 1, 8), filename, ".beams"))
		}
	# Define the output:
	bgns = c(bgns, nrnp, nrna, hins, list(fltr=filter, nsdn=nsdn, sfsi=surface))
	# Write the noise estimates:
	if(!is.null(con)){
		if(!isTRUE(file.info(dirname(con))$isdir)){
			dir.create(dirname(con), recursive=TRUE)
			}
		write.TSD(bgns, con=con, dimension=TRUE, ...)
		}
	# Output the noise estimates:
	bgns
	##################################################
	##################################################
	}
