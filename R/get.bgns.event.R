#*********************************************
#*********************************************
#' Estimates the mean of exponentially distributed background noise for each beam of an underwater acoustic system, based on the lower fraction of the values specified by 'prob', and writes the result to file in the TSD format if 'con' is given.
#'
#' @param con  is the connection object or a character string naming the output file.
#' @param event  is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
#' @param cruise  is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param t  is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
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
#' @param pdns_scale  is used for scaling the data to ensure best performance of the optimization, and corersponds to the mean of the data (default if set to NULL).
#' @param lowerpar,startpar,upperpar are the lower, start and upper parameter values representing the magnitude, kurtosis and phase, respectively, where the magnitude is given relative to 'pdns_scale'.
#' @param kern  is either a vector of wieghts summing to 1, specifying the smoothing of the spectral frequency used to identify the harmonics, or a single numeric = 1, specitying a Gaussian kernel with sd=1.
#' @param high  is used in resampleHigh() in get.bgns.event() to reduce the impact from unusually high values in the fft.
#' @param ...  variables used in write.TSD() (such as 'ow' for overwriting existing files).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD arr.ind2ind ind.expand NAs read.TSDs write.TSD zeros
#' @importFrom stats dnorm fft filter optim rnorm runif
#'
#' @export
#' @rdname get.bgns.event
#'
get.bgns.event<-function(con=NULL, event=NULL, cruise=2009116, t=NULL, esnm="MS70", treat  = c("bgns", "pdns"), noise=NULL, fftthreshold=0.145, nsind=0.75, acfq=100, prob=0.2, type=6, tries=c(10,10), thr=1e-3, trim=0.1, search=4, kurtosis=c(0.1,3), reffan=16, pdns_scale=1e-14, lowerpar=c(1e-3,0.001,-2*pi), startpar=c(1,2,3), upperpar=c(1e3,6,4*pi), kern=NULL, high=Inf, ...){
	
	############### LOG: ###############
	# Start: 2010-03-12 - Clean version.
	# Update: 2012-05-23 - Changed to estimate the periodic noise as well, using optim() and estimating the background noise and periodic noise simultaneously.
	# Update: 2012-05-23 - Changed to estimate the background noise first, then optimize to estimate the periodic noise.
	# Update: 2012-05-23 - Changed method for estimation of the background noise to using meanexp.quantile() instead of meanexp[dot]conditional().
	# Update: 2012-05-29 - Expanded to several events.
	# Update: 2012-06-13 - Changed to identify the beams affected by the periodic noise by storing the relative height of the highest peak in spectral density compared to the specral density of the zero-frequency, averaging over all pings, and thresholding with 'fftthreshold'.
	# Update: 2012-07-06 - Changed to 10^() instead of exp().
	# Update: 2012-07-13 - Radical change in the estimation of the background noise and periodic noise for the periodic beams. Now these are estimated iteratively by alternating between estimation of the background noise and the periodic noise and evaluating the expression Q(u) = 1 - sum(exp(-u/(beta_B + pulsenoise)). Also added the paramtters tries=10 and thr=1e-3, specifying the number of tries in the iteration, and the threshold for the absoulte difference between the previous and the current estimate, divided by the current estimate, at which the iteration stops.
	# Update: 2012-09-21 - Changed the estimation of the harmonics 'harm' and the phase 'ps3M' so that these are constant within each fan, according to the assumptions made in Papaer II in the PhD.
	# Update: 2012-09-23 - Restructured the iteration and added comments. Also the function values are stored and used when the number of iterations exceeds the value of 'search'. Expanded the input 'tries' to be of length 2, specifying the number of tries for the iteration [1] and the optimization [2].
	# Update: 2012-10-26 - Removed defaults for 't', 'event' and 'cruise'.
	# Update: 2012-10-29 - Added the option 'lowerkurtosis', which discards the estimations for which the kurtosis estimate is lower than 0.5, to avoid erroneous estimates from affecting the mean estimates.
	# Update: 2012-10-31 - Added the options 'pdns_scale', 'lowerpar', 'startpar' and 'upperpar'.
	# Update: 2012-11-08 - After extensive testing the following issues were raised: (1) The harmonics estimated from the frequency spectrum is sensitive to the length of the beam segments. This problem may be reduced solved by smoothing the frequency spectrum by a gaussian kernel (sd=1 for n=1000), which is added as an option ('kern') in the new version. (2) The periodc noise in beams 351:375 and 401:425 differ from the ordinary periodic noise by appearing to have two peaks, one on either side of whre the peak is expected. The kurtosis is restricted to be >lowerkurtosis to avoid not finding the periodic noise, and for these beams one of the peaks in the pairs of peaks is typically selected. (3) An error in the estimation of the mean phase parameters was detected, where the weighting should be pns1 * 10^pns2. (4) The 'fftthreshold' was set to 0.1 if the kernel is given as a single numeric = 1, indicating a Gaussian kernel with sd=1. This was found from multiple testing on 10 consecutive time steps for different values of 'fftthreshold': nperiodic=cbind(seq(0.005,0.14,0.005),c(rep(500,12),425,95,68,68,68,67,66,66,66,63,60,59,59,55,54,54)); ploto(nperiodic). However, the value was tuned to the nice value 0.1 after processing all the 100 time steps used in the noise paper.
	# Update: 2012-11-19 - Changed the name of 'lowerkurtosis' to 'kurtosis', and added the option of both lower and upper limit for the kurtosis in this parameter.
	# Update: 2012-11-24 - Added the simple mean of each beam, regardless of periodic noise, to the output.
	# Last: 2013-07-19 - Removed 'near.range' and renamed 'ind1' to 'nsind'.
	
	
	# Return an empty list if any of 't' or 'event' are empty:
	if(length(unlist(t))==0 || length(unlist(event))==0){
		return()
	}
	
	# Transform 't' to list:
	if(length(t)>0 && !is.list(t)){
		t = list(t)
	}
	nevents = length(event)
	
	# Check for compatibility between events if more than one event is given:
	if(length(t) != nevents){
		stop("Unequal number of events and time vectors given for the backgound noise")
	}
	if(length(cruise) != nevents){
		cruise = rep(cruise, length.out=nevents)
	}
	
	# Read the length of the sample intervals, the number of beams, and the length of the beams:
	sint = NULL
	numb = NULL
	lenb = NULL
	# Read the time steps of the events:
	indt = vector("list", nevents)
	dimb = vector("list", nevents)
	utim = NULL
	
	# Run through the events:
	for(ev in seq_len(nevents)){
		data = read.event(event=event[ev], cruise=cruise[ev], var="beams", esnm=esnm)
		sint = c(sint, data$sint)
		numb = c(numb, data$numb)
		lenb = c(lenb, max(data$lenb))
		time = read.event(event=event[ev], cruise=cruise[ev], t=t[[ev]], var="time", esnm=esnm)
		indt[[ev]] = time$indt
		dimb[[ev]] = c(data$numb/length(unique(data$freq)), length(unique(data$freq)))
		utim = c(utim, time$utim)
	}
	
	# Define the number of time steps:
	numt = length(unlist(indt))
	if(numt==0){
		stop("No valid time steps selected")
	}
	
	# Check whether the beam configuration is equal across events:
	if(!all(sint==sint[1]) || !all(numb==numb[1]) || !all(lenb==lenb[1])){
		stop("Events must have the same beam configuration")
	}
	else{
		sint = sint[1]
		numb = numb[1]
		lenb = lenb[1]
		dimb = dimb[[1]]
	}
	
	# Interpret the valid voxels 'nsind' (see description of variables):
	nsind = ind.expand(nsind, lenb, drop=TRUE)
	
	
	
	
	
	
	
	
	if("pdns" %in% treat) {
		# Mode function:
		getharm = function(x){
			x = as.numeric(names(which.max(table(x))))
			if(length(x)==0){
				1
			}
			else{
				x
			}
		}
		# Function for calculating the mean of angles:
		angleMean = function(ang, w=1){
			atan2(sum(w*sin(ang), na.rm=TRUE), sum(w*cos(ang), na.rm=TRUE)) %% (2*pi)
		}
		
		if(identical(kern, 0)){
			kern = NULL
		}
		if(length(kern)==1){
			if(kern==1){
				fftthreshold = 0.1 # Gives slightly fewer periodic beams than the unfiltered spectral density with fftthreshold=0.145
			}
			kern = dnorm(seq(-3*kern,3*kern), 0, kern)
			#kern = dnorm(-3:3,0,1)
			kern = kern/sum(kern)
		}
		tries = rep(tries, length.out=2)
		
		
		
		
		
		# Check if the reference fan index is valid for all events:
		for(ev in seq_len(nevents)){
			if(!reffan%in%seq_len(dimb[[ev]][2])){
				stop(paste("reference fan", reffan, "chosen to high"))
			}
		}
		
		
		
		# Define the factor in the sine wave (sin(ax) gives period dt = 2*pi/a => a = 2*pi/dt. dt = time between peaks 1/acfq divided by pulselength 'sint', and divided by 2 to fit to the observed frequencies => a = 2*pi * acfq*sint):
		twopif = 2*pi * acfq*sint
		
		# Define the function that calculates the periodic noise added the background noise:
		pulsenoise = function(par, t, n, twopif, bgns){
			bgns + par[1]*10^(par[2]*sin(n*twopif*t + par[3]))
		}
		
		# Define the function used in the opimazition (bgns + pns1 * exp(pns2 * sin(n*a*j - pns3))):
		pulsenoiseSum = function(par, x, t, n, twopif, bgns){
			#sum(abs(pulsenoise(par,t,n,twopif,bgns)-x)) No longer in use
			sum((pulsenoise(par, t, n, twopif, bgns)-x)^2)
		}
		
		# Define the difference between the lower and the upper parameter values:
		diffpar = upperpar-lowerpar
		
		
		
		##### Execution and output #####
		tol = 1e-5
		bg0M = NAs(numb, numt)
		bgnM = NAs(numb, numt)
		pn1M = NAs(numb, numt)
		pn2M = NAs(numb, numt)
		pn3M = NAs(numb, numt)
		fnvb = NAs(numb, numt)
		fnvp = NAs(numb, numt)
		fnvt = NAs(numb, numt)
		ntry = NAs(numb, numt)
		NTRY = NAs(numb, numt)
		
		if(length(noise)>0 && is.character(noise)){
			noise = read.TSDs(noise, var=c("harm","badb"))
		}
		# Get the harmonics and indices for the bad beams, and if these are missing in 'noise' it results in re-generating the variables:
		suppressWarnings(harm<-rep(noise$harm, length.out=numb))
		badb = rep(as.logical(noise$badb), length.out=numb)
		
		
		# If the harmonics are not present, generate from the data:
		if(length(harm)==0){
			
			harM = NAs(numb, numt)
			secondmax = NAs(numb, numt)
			
			### 1. Apply the Fast Fouriertransform (fft) to the beams, and identify the beams that are affected by the periodic noise: ###
			secondmax = NAs(numb, numt)
			p = 1
			
			
			for(ev in seq_len(nevents)){
				cat(paste("Checking for periodic noise, event ", event[ev], "\nPings:\n"))
				
				for(i in indt[[ev]]){
					cat(i, " ", sep="")
					# Import acoustic data and ignore near-range and apply the subset specified by 'nsind':
					Sv = read.event(event=event[ev], cruise=cruise[ev], t=i, var="vbsc", TVG=FALSE, esnm=esnm)$vbsc
					thisvalidvoxels = intersect(which(!apply(Sv, 1, function(x) any(is.na(x)))), nsind)
					if(length(pdns_scale)==0 && ev==1 && i==indt[[ev]][1]){
						pdns_scale = mean(Sv[thisvalidvoxels,])
					}
					
					Sv = Sv[thisvalidvoxels,,drop=FALSE]/pdns_scale
					
					# Identify the periodic noise:
					for(j in seq_len(numb)){
						fs = 1/(sint*acfq)
						# Define the values that will not be used when searching for peaks in the frequency spectrum:
						discard = seq_len(1/(fs/length(thisvalidvoxels))/2)
						ndiscard = length(discard)
						# Resample the values which are higher that 'high', as divided by the trimmed mean, then Fourier transform:
						fftBeam = fft(resampleHigh(Sv[,j], high, trim=0.1))
						m = Mod(fftBeam)[seq_len(length(fftBeam)/2)]
						maxm = max(m)
						# Smooth the fourier spectrum:
						if(length(kern)>0){
							m = filter(m, kern)
						}
						secondmax[j,p] = max(m[-discard], na.rm=TRUE)/maxm
						# Remove the first element, which represents the infinite period (mean):
						harM[j,p] = (which.max(m[-discard])+ndiscard)*fs/length(thisvalidvoxels)
					}
					p = p+1
				}
				cat("\n")
			}
			
			# Identify the periodic beams by those who have mean relative sprectral density for the highest peak above the threshold 'fftthreshold':
			rsdM = secondmax
			rspd = rowMeans(secondmax)
			badb = rspd>fftthreshold
			
			# Round off to the nearest integer 1 or 2:
			hM12 = harM
			hM12[hM12<=0.5] = NA
			hM12[hM12<=1.5 & hM12>0.5] = 1
			hM12[hM12<=2.5 & hM12>1.5] = 2
			hM12[hM12>2.5] = 1
			
			# Set all values of the non-periodic beams to NA:
			hM12[badb==0,] = NA
			
			# Average and round off for each beam:
			dim(hM12) = c(dimb, numt)
			harm = apply(hM12, 2, getharm)
			harm[colSums(array(badb, dimb))==0] = 1
			dim(hM12) = c(prod(dimb), numt)
			# Expand the estimates of 'harm' into the full dimensions of the sonar:
			harm = rep(harm, each=dimb[1])
		}
		else{
			rsdM = NA
			harM = NA
			hM12 = NA
			rspd = NA
			if(length(badb)==0){
				badb = harm>0
			}
		}
		
		
		if(sum(badb)>0){
			indbadbeams = which(badb)
		}
		else{
			indbadbeams = NULL
		}
		
		
		### 2. Estimate the background noise 'beta_N' by the simple mean for the healthy beams and by the mean of the lower quantile of the data for beams affected by the periodic noise: ###
		p = 1
		for(ev in seq_len(nevents)){
			cat(paste("Estimating background noise, event ", event[ev], "\nPings:\n", sep=""))
			
			for(i in indt[[ev]]){
				cat(i," ",sep="")
				# Import acoustic data and ignore near-range and apply the subset specified by 'nsind':
				Sv = read.event(event=event[ev], cruise=cruise[ev], t=i, var="vbsc", TVG=FALSE, esnm=esnm)$vbsc
				thisvalidvoxels = intersect(which(!apply(Sv, 1, function(x) any(is.na(x)))), nsind)
				Sv = Sv[thisvalidvoxels,,drop=FALSE]/pdns_scale
				
				# (a) For the healthy beams use the simple mean:
				bg0M[,p] = colMeans(Sv)
				#bg0M[,p] = meanexp.quantile(Sv, prob = prob, type = type, MARGIN = 2)
				bgnM[!badb,p] = bg0M[!badb,p]
				
				# (b) Calculate the parameter of the exponential distribution by use of the lower portion of the values:
				# Estimate the expectation by the quantile divided by ln(1/(1-p)), derived from F(x)=p:
				# Previously meanexp[dot]conditional() was used, but this included a bias.
				if(length(indbadbeams)>0){
					initialbgnM = meanexp.quantile(Sv[,indbadbeams,drop=FALSE], prob=prob, type=type)
					
					# Optimize each of the beams identified to be affected by the periodic noise:
					for(j in seq_along(indbadbeams)){
						# Define the difference between the previous and the current background noise estimate, used for ternimating the iteration between estimation of periodic noise and background noise:
						diff = Inf
						# Define the number of iterations:
						TR = 0
						
						# Define objects storing the estimates of the background noise and the periodic noise, and objects storing the function values at the optimum values:
						thisbgns = initialbgnM[j]
						thispns1 = NA
						thispns2 = NA
						thispns3 = NA
						valuebgns = NA
						valuepdns = NA
						sumofsquares = NA
						
						# Loop through the iterations and terminate if the background noise estimate does not change more than 0.001 of the curent value, or if the number of iterations has exceeded the limit:
						while(diff>thr  && TR<tries[1]){
							TR = TR+1
							# Optimize the perioic noise and append the estimates to the vectors of the estimates for the current beam of the current time step:
							o = optim(startpar, pulsenoiseSum, x=Sv[,indbadbeams[j]], t=thisvalidvoxels, n=harm[indbadbeams[j]], twopif=twopif, bgns=thisbgns[TR], method="L-BFGS-B", lower=lowerpar, upper=upperpar)
							thispns1 = c(thispns1, o$par[1])
							thispns2 = c(thispns2, o$par[2])
							thispns3 = c(thispns3, o$par[3])
							
							# Check if the value of the kurtosis is too close to the limits:
							tr = 0
							
							while(any((thispns2[TR+1]-lowerpar[2])<tol*diffpar[2], (upperpar[2]-thispns2[TR+1])<tol*diffpar[2]) && tr<tries[2]){
								tr = tr+1
								# Draw random starting parameters:
								uu = c(startpar[1], runif(1)*diffpar[2]+lowerpar[2], startpar[3])
								# Re-optimize and store the estimates at the current position:
								o = optim(uu, pulsenoiseSum, x=Sv[, indbadbeams[j]], t=thisvalidvoxels, n=harm[indbadbeams[j]], twopif=twopif, bgns=thisbgns[TR], method="L-BFGS-B", lower=lowerpar, upper=upperpar)
								thispns1[TR+1] = o$par[1]
								thispns2[TR+1] = o$par[2]
								thispns3[TR+1] = o$par[3] %% (2*pi)
							}
							
							# Store the number of tries and the final function value of the optimization:
							ntry[indbadbeams[j],p] = tr
							valuepdns = c(valuepdns,o$value)
							
							# Calculate the periodic noise for use in the estimation of the background noise:
							a = pulsenoise(c(thispns1[TR+1], thispns2[TR+1], thispns3[TR+1]), t=thisvalidvoxels, n=harm[indbadbeams[j]], twopif=twopif, bgns=0)
							
							# Re-estimate the background noise and store the estimate and the function value:
							m = meanexp.quantile_a(Sv[,indbadbeams[j],drop=FALSE], lower=0, upper=100, prob=prob, a=a, type=type)
							valuebgns = c(valuebgns, m$objective)
							thisbgns = c(thisbgns, m$minimum)
							# Store the sum of squares of the estimated noise:
							sumofsquares = c(sumofsquares, sum((Sv[,indbadbeams[j],drop=FALSE]-m$minimum-a)^2))
							
							# Update the difference value:
							diff = abs(thisbgns[TR+1]-thisbgns[TR])/thisbgns[TR+1]
						}
						NTRY[indbadbeams[j],p] = TR
						
						# If the number of tries exceed 5, extract the value with the lowest sum of the rank of the function values for both estimations:
						if(TR>=search){
							#thisrank = rank(valuebgns)+rank(valuepdns)
							#thisrank = max(which(thisrank==min(thisrank)))
							thisrank = which.min(sumofsquares)
						}
						else{
							thisrank = TR+1
						}
						# Insert the estimates to the output:
						bgnM[indbadbeams[j],p] = thisbgns[thisrank]
						pn1M[indbadbeams[j],p] = thispns1[thisrank]
						pn2M[indbadbeams[j],p] = thispns2[thisrank]
						pn3M[indbadbeams[j],p] = thispns3[thisrank]
						fnvb[indbadbeams[j],p] = valuebgns[thisrank]
						fnvp[indbadbeams[j],p] = valuepdns[thisrank]
						fnvt[indbadbeams[j],p] = sumofsquares[thisrank]
					}
				}
				p = p+1
			}
			cat("\n")
		}
		
		bgnM = bgnM*pdns_scale
		bg0M = bg0M*pdns_scale
		pn1M = pn1M*pdns_scale
		
		# Preapare variables for writing:
		validpings = kurtosis[1]<pn2M & pn2M<kurtosis[2]
		validpings[is.na(validpings)] = FALSE
		if(numt==1){
			dim(validpings) = c(numb, 1)
		}
		
		# Estimate the mean magnitude and mean kurtosis of the periodic noise, and the mean background noise, by discarding the estimates for which the kurtosis is below the kurtosis interval 'kurtosis':
		pns1 = NAs(numb)
		pns2 = NAs(numb)
		# Average over the time steps:
		bgns = rowMeans(bgnM)
		bgn0 = rowMeans(bg0M)
		for(i in which(badb)){
			pns1[i] = mean(pn1M[i,validpings[i,]], trim=trim, na.rm=TRUE)
			pns2[i] = mean(pn2M[i,validpings[i,]], trim=trim, na.rm=TRUE)
			#bgns[i] = mean(bgnM[i,validpings[i,]],trim=trim,na.rm=TRUE)
		}
		
		# Equal the phase estimates for each fan:
		pn30 = pn3M
		# Reshape into [I1,I2,P]
		dim(pn1M) = c(dimb, numt)
		dim(pn2M) = c(dimb, numt)
		dim(pn3M) = c(dimb, numt)
		dim(validpings) = c(dimb, numt)
		# For each fan of each time step, calculate the mean angle of the phase angles, using the magnitude estimates as the lengths of vectors in the direction of the phase angles:
		for(i2 in seq_len(dim(pn3M)[2])){
			for(p in seq_len(dim(pn3M)[3])){
				pn3M[,i2,p] = angleMean(pn3M[,i2,p][validpings[,i2,p]], pn1M[,i2,p][validpings[,i2,p]] * 10^(pn2M[,i2,p][validpings[,i2,p]]))
			}
		}
		# Also, for each fan, calculate the mean angle of the phase angles, using the magnitude estimates as the lengths of vectors in the direction of the phase angles:
		pns3 = NAs(dimb)
		for(i2 in seq_len(dim(pn3M)[2])){
			pns3[,i2] = angleMean(pn3M[,i2,][validpings[,i2,]], pn1M[,i2,][validpings[,i2,]] * 10^(pn2M[,i2,][validpings[,i2,]]))
		}
		dim(pns3) = NULL
		
		# Reset dimensions:
		dim(pn1M) = c(prod(dimb), numt)
		dim(pn2M) = c(prod(dimb), numt)
		dim(pn3M) = c(prod(dimb), numt)
		
		# Since the correlation depends on where the peaks occur, we convert to the locations of the peaks of the periodic noise for each fan, using the most dominant fan (16 for the MS70) as a reference, and select the location of the first peak preceding the reference peak, of each of the other fans. Converting for locations of the peaks is done by solving 2*pi*f_P * j + phase = z*(2*pi) + pi/2, for j, which results in j = (z*(2*pi) + pi/2 - phase) / (2*pi*deltaT*f_AC * n), where f_P = 2*pi*detlaT*f_AC is the frequency of the periodic noise, deltaT is the sample interval duration, f_AC is the frequency of the alternating current, and j is the sampling interval index.
		
		# Define a matrix of one row for each fan with time steps along the rows, holding the locations of the peaks of the periodic noise:
		peaksoffans = zeros(dimb[2], numt)
		# Define the difference from the peak values of a fan to the peak values of the same time steps for the reference fan:
		difftoreffan = zeros(dimb[2], numt)
		# Define the mean difference in peak location from a fan to the reference fan:
		meandiff = zeros(dimb[2])
		# Define the new estimates of the phase angles calculated by the mean difference in peak location to the reference fan:
		pn3I = zeros(dimb[2], numt)
		
		# Define a matrix of indexes identifying the fans:
		fans = array(seq_len(numb), dim=dimb)	
		# For the reference fan calculate the location of the peaks for four different values of 'z':
		thispeaks = outer(pi/2 - pn3M[fans[1,reffan],], 0:3 *2*pi,"+") / (2*pi*sint*acfq * harm[fans[1,reffan]])
		# Extract the first positive peak location for each time step:
		peaksoffans[reffan,] = apply(thispeaks, 1, function(x) min(x[x>0]))
		
		# For the other periodic fans (12, ..., 15 and 17) calculate the peak location closest to the peaks in the reference fan:
		badb = as.numeric(badb)
		otherPeriodicFans = which(apply(array(badb, dim=dimb), 2, function(x) sum(x)>0))
		for(i in otherPeriodicFans){
			# For the current fan calculate the location of the peaks for 7 different values of 'z' (probably more than enough):
			thispeaks = outer(pi/2 - pn3M[fans[1,i],], 0:6 *2*pi,"+") / (2*pi*sint*acfq * harm[fans[1,i]])
			# Get the difference between the peak locations of the current fan and the reference fan:
			absdiff = abs(thispeaks-peaksoffans[reffan,])
			# Extract the position of the closest value:
			absdiff = arr.ind2ind(cbind(seq_len(nrow(thispeaks)), apply(absdiff,1,which.min)),dim(thispeaks))
			# Extract the peak location closest to the peak locations of the reference fan:
			peaksoffans[i,] = thispeaks[absdiff]
			# Calculate the difference in the peak location between the current fan and the reference fan
			difftoreffan[i,] = peaksoffans[i,]-peaksoffans[reffan,]
			
			# Transpose the peak locations that are close to one period fathrer than the closest locations:
			w = 1/(sint*acfq*harm[fans[1,i]])
			if(diff(range(difftoreffan[i,])) > (w-0.5)){
				cat(paste("Transposing value in fan ",i," to prepare for extracting the mean difference between peaks","\n",sep=""))
				mid = max(difftoreffan[i,])-w/2
				difftoreffan[i,][difftoreffan[i,]>mid] = difftoreffan[i,][difftoreffan[i,]>mid] - w
			}
			
			# Calculate the mean difference in peak location to the reference fan:	
			meandiff[i] = mean(difftoreffan[i,])
			
			# Convert back to phase angle for each time step. This is obtained by solving 2*pi*f_P * j + phase = z*(2*pi) + pi/2, for phase, which results in phase = z*(2*pi) + pi/2 - phase - 2*pi*deltaT*f_AC * n * j, where f_P = 2*pi*detlaT*f_AC is the frequency of the periodic noise, deltaT is the sample interval duration, f_AC is the frequency of the alternating current, and j is the sampling interval index.:
			pn3I[i,] = (pi/2-2*pi*sint*acfq * harm[fans[1,i]] * (peaksoffans[reffan,]+meandiff[i])) %% (2*pi)
			
			# Convert back to phase angle averaged over time steps:
			pns3[fans[,i]] = (pi/2-2*pi*sint*acfq * harm[fans[1,i]] * meandiff[i]) %% (2*pi)
		}
		
		# Expand into desired dimensions:
		pn3I = rep(pn3I,each=dimb[1])
		dim(pn3I) = c(prod(dimb),numt)
		
		# Define output values for the difference in peak location between all beams and the reference fan of all time steps:
		DPrf = rep(difftoreffan,each=dimb[1])
		dim(DPrf) = c(prod(dimb),numt)
		# Define output values for the mean over all time steps of difference in peak location between all beams and the reference fan:
		mDrf = rep(meandiff,each=dimb[1])
		dim(mDrf) = c(prod(dimb))
	}
	
	
	
	
	
	
	
	
	
	
	
	
	else{
		
		bg0M = NAs(numb, numt)
		
		p = 1
		for(ev in seq_len(nevents)){
			cat(paste("Estimating background noise, event ", event[ev], "\nPings:\n", sep=""))
			
			for(i in indt[[ev]]){
				cat(i," ",sep="")
				# Import acoustic data and ignore near-range and apply the subset specified by 'nsind':
				Sv = read.event(event=event[ev], cruise=cruise[ev], t=i, var="vbsc", TVG=FALSE, esnm=esnm)$vbsc
				#thisvalidvoxels = intersect(which(!apply(Sv, 1, function(x) any(is.na(x)))), nsind)
				thisvalidvoxels <- 301:1310
				Sv = Sv[thisvalidvoxels,,drop=FALSE]/pdns_scale
				
				# (a) For the healthy beams use the simple mean:
				bg0M[,p] = colMeans(Sv)
				#bg0M[,p] = meanexp.quantile(Sv, prob = prob, type = type, MARGIN = 2)
				
				p = p+1
				}
				
			}
			cat("\n")
		
		
		bg0M = bg0M*pdns_scale
		bgnM <- bg0M
		# Average over the time steps:
		bgns = rowMeans(bgnM)
		bgn0 = rowMeans(bg0M)
		
		
		badb = NAs(numb)
		pns1 = NAs(numb)
		pns2 = NAs(numb)
		pns3 = NAs(numb)
		harm = NAs(numb)
		rspd = NAs(numb)
		rsdM = NAs(numb, numt)
		pn1M = NAs(numb, numt)
		pn2M = NAs(numb, numt)
		pn3M = NAs(numb, numt)
		pn3I = NAs(numb, numt)
		pn30 = NAs(numb, numt)
		harM = NAs(numb, numt)
		hM12 = NAs(numb, numt)
		DPrf = NAs(numb, numt)
		mDrf = NAs(numb)
		#acfq = acfq, 
		LOWP = lowerpar
		BEGP = startpar
		UPPP = upperpar
		SCLE = pdns_scale
		fnvb = NAs(numb, numt) 
		fnvp = NAs(numb, numt) 
		fnvt = NAs(numb, numt) 
		ntry = NAs(numb, numt) 
		NTRY = NAs(numb, numt) 
		FTth = fftthreshold
		pn2I = kurtosis
		
		#utim = utim, 
		prex = prob
		#esnm = esnm
	}
	
	
	# Write and return the data:
	out = list(bgns=bgns, badb=badb, pns1=pns1, pns2=pns2, pns3=pns3, harm=harm, rspd=rspd, rsdM=rsdM, bgnM=bgnM, bgn0=bgn0, bg0M=bg0M, pn1M=pn1M, pn2M=pn2M, pn3M=pn3M, pn3I=pn3I, pn30=pn30, harM=harM, hM12=hM12, DPrf=DPrf, mDrf=mDrf, utim=utim, acfq=acfq, LOWP=lowerpar, BEGP=startpar, UPPP=upperpar, SCLE=pdns_scale, fnvb=fnvb, fnvp=fnvp, fnvt=fnvt, ntry=ntry, NTRY=NTRY, FTth=fftthreshold, prex=prob, pn2I=kurtosis, esnm=esnm)
	if(!is.null(con)){
		write.TSD(out,con=con,numt=numt,...)
		}
	out
}
