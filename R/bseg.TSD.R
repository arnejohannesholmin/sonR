#*********************************************
#*********************************************
#' Estimates the probability of not school in voxels of sonar data.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param t  is either the indexes of the pings to be treated, as listed from 1 to the number of pings in the event, or the time point given as a string "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF".
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event.
#' @param bgns  is an optional list of noise estimates (must contain the background noise 'bgns', and may contain the periodic noise estimates 'pns1', 'pns2' and 'harm'). The phase parameter 'pns2' is extracted from the data. Use 'bgns' to suppress adding periodic noise or near range noise in the noise which is compared to the data, by setting bgns=list(acfq=NULL) or bgns=list(nr0a=NULL).
#' @param beta0  is the minimum schooling threshold for the volume backscattering coefficient. For voxels where beta0<noise: beta0=noise.
#' @param beta1  is the maximum schooling threshold for the volume backscattering coefficient, defining the probability distribution of the signal. Should be chosen on the basis of the maximum packing density of the observed species.
#' @param factor  is one or several factors to multiply the schooling thresholds 'beta0' and 'beta1' by, increasing the memory occupied by the function by 50 \% of the original memory (factor==NULL) for each element of 'factor'. This option was included to reduce CPU time in the case that both an unbiased and an enlarged segmentation mask is requested.
#' @param ind  is a vector or list of indexes along the beams, as input to ind.expand(), used to select the subset over which the estimation of background noise is done. Also voxels not included in this subset are assigned p=1.
#' @param nsind  is a vector of indexes along the beams, as input to ind.expand(), used to select the subset over which the estimation of the phase of the periodic noise is done. If given as a single numeric, the outermost 'nsind' voxels are used in each beam.
#' @param smind  is a list of indexes used to define the subset which is smoothed spatially to accumulate the probabilities. The default avoids smoothing densely packed voxels, which reduces computational time.
#' @param h  is the bandwidth of the spatial Gaussian kernel smoother, given as the standard deviation of the Gaussian distribution. For no smoothing use h=NULL.
#' @param alpha  is the significance level of the hypothesis testing of H0: school not present in the voxel, against H1: school present in the voxel. To return non-thresholded data, use alpha=NULL.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param hins_add  is the number of voxels that should be discarded on both sides of high intensity noise voxels voxels along beams, used for accounting for possible high values that are related to the high intensity noise but not classified as such voxels.
#' @param phase  is FALSE if any of 'pn3M' (phase for each time step) or 'pns3' (phase equal for all time steps) given in 'bgns' or read from the noise file located by the funciton noise.path.event() should be used, as oposed to estimating the phase from the data for each time step. This is only recommended for simulated data where the phase is constant over all time steps, and saves some CPU time.
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param subtractNoise  is TRUE if the original acoustic data should be subtracted noise, which is used when returning the total and mean volume backscatter from the segment.
#' @param adds  is an optional list of variables overriding the variables in read from the event.
#' @param sim  is a TRUE if smoothing should be done only along the first dimensions, simultaneously over the stages of the last dimension. If 'sim' is an integer larger than 1, the positions 'coords' are used 'sim' times, and the data 'x' should have length 'sim' times the length of one coordinate of 'coords'.
#' @param na.rm  is single integer representing the dimension along which NAs are discarded from the smoothing in the case that sim==1. For example, if na.rm=2 and the dimension of 'x' is [5,12,7], and x[3,2:4,5]=NA, then all data x[,2:4,] will be excluded from the smoothing and set to NA. If na.rm=FALSE, no NAs should be contained in the data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom expint expint_E1
#' @importFrom SimradRaw apply.TVG
#' @importFrom TSD arr.ind2ind ind.expand NAs ones read.TSDs zeros
#'
#' @export
#' @rdname echoIBM.vbsc2p.event
#'
echoIBM.vbsc2p.event <- function(event=1, t=1, cruise=2009116, bgns=NULL, beta0=NULL, beta1=NULL, factor=NULL, ind=list(-(1:30),NULL), nsind=0.75, smind=list(-(1:300)), range=list(), subset=NULL, h=NULL, alpha=NULL, dir.data=NULL, hins_add=10, phase=TRUE, TVG.exp=2, esnm="MS70", subtractNoise=TRUE, adds=list(), sim=TRUE, na.rm=1, allow.old=FALSE, TOV=0){
	
	##### Preparation #####
	# A small funciton used to get the minimum of each row in a two column matrix (much faster than apply(x,1,min)):
	rowMin2 <- function(x){
		a <- x[,2] > x[,1]
		x[a,2] <- x[a,1]
		x[,2]
	}
	
	# Define the function for calculating the cummulative distribution. This function takes up to one second for one ping of MS70 data, mostly due to the expint_E1() function:
	P <- function(beta0, beta1, x, betaN){
		
		# 2024-03-14: THIS was not correct. Ann exponential integral should be used:
		# # The integrand which is the exponential CDF divided by beta1, which is the upper limit of the uniform prior:
		# integrand <- function(b, beta1, x, betaN){
		# 	1/beta1 * exp(-x / (b + betaN))
		# }
		# 
		# H <- function(beta0, beta1, x, betaN){
		# 	integrand(b=beta0, beta1=beta1, x=x, betaN=betaN) - integrand(b=0, beta1=beta1, x=x, betaN=betaN)
		# }
		# 
		# out <- H(beta0) / H(beta1)
		# 
		# out[outsidemaxarg_0] <- .Machine$double.eps
		# out
		
		H <- function(b, x, betaN, beta1){
			expint::expint_E1(x / (b + betan)) - expint::expint_E1(x / ( betan))
		}
		
		H(beta0) / H(beta1)
		
	}
	
	
	lf <- length(factor)
	numt <- length(t)
	
	# 'beta0' or 'beta1' should be of the same length as the number of time steps:
	if(!length(beta0)) {
		stop("beta0 must be given, either as a single value or a vector (repeated to the number of time steps)")
	}
	if(!length(beta1)) {
		stop("beta1 must be given, either as a single value or a vector (repeated to the number of time steps)")
	}
	if(length(beta0)!=numt){
		beta0 <- rep(beta0,length.out=numt)
	}
	if(length(beta1)!=numt){
		beta1 <- rep(beta1,length.out=numt)
	}
	
	# Update 'esnm':
	esnm <- read.event(event=event, cruise=cruise, esnm=esnm, var="esnm", msg=FALSE)$esnm
	
	# 'sim' can only be TRUE for sonars:
	if(!is.sonar(esnm[[1]])){
		sim <- FALSE
	}
	
	# For simultaneous smoothing for all time steps using the voxels positions in the coordinate system of the vessel, read the voxels positions only once:
	#if(sim){
	#	dynamicvar=c("vbsc","time")
	#	staticvar=c("psxx","psyx","pszx","volx","beams")
	#	cs="v"
	#}
	#else{
	#	dynamicvar=c("vbsc","time","esnm","psxx","psyx","pszx")
	#	staticvar=c("volx","beams")
	#	cs="g"
	#}
	if(sim){
		dynamicvar=c("vbsc", "time", "beams")
		staticvar=c("psxx", "psyx", "pszx", "volx")
		cs="v"
	}
	else{
		dynamicvar=c("vbsc", "time", "esnm", "psxx", "psyx", "pszx", "beams")
		staticvar=c("volx")
		cs="g"
	}
	noisevar=c("hini","hins")
	
	#browser()
	
	# Read the acoustic data:
	data <- read.event(event=event, cruise=cruise, esnm=esnm, t=t, var=dynamicvar, msg=FALSE, cs=cs, allow.old=allow.old, TOV=TOV)
	# Add high intensity noise:
	temp <- read.event(event=event, cruise=cruise, t=t, var=noisevar, msg=FALSE, allow.old=allow.old)
	data[names(temp)] <- temp
	# Read the voxel position data in the coordinate system of the vessel, which is static for all time steps, allowing for the smoothing to be done over several time steps at the time. Also add the sonar configuration:
	temp <- read.event(event=event, cruise=cruise, t=1, var=staticvar, msg=FALSE, cs=cs, allow.old=allow.old)
	data[names(temp)] <- temp
	# Store the original dimension of the acoustic data:
	if(length(dim(data$vbsc))==2){
		dim(data$vbsc) <- c(dim(data$vbsc),numt)
	}
	olddim <- dim(data$vbsc)
	# Read background and periodic noise data, added to the optional list 'bgns':
	noisefiles <- noise.path(cruise=cruise,event=event,esnm=esnm,dir.data=dir.data)
	data <- c(adds, data, bgns, read.TSDs(noisefiles,var=c("bgns","badb","pns1","pns2","pns3","pn3M","harm","acfq","nr0a")))
	periodic <- which(data$badb==1)
	rm(bgns)
	gc()
		
	# If zeros occur in data$bgns, add a small value:
	if(min(data$bgns)==0){
		data$bgns[data$bgns==0] <- .Machine$double.eps
	}
	
	# If the required information is present in 'data', estimate the phase of the periodic noise at each time step, from a portion of the sonar volume dominated by noise, or alternatively, not occupied by schools:
	if(phase){
		data$pn3M <- get.pdns_phase.event(event=event, cruise=cruise, esnm=esnm, t=t, noise=data, nsind=nsind)$pn3M
	}
	
	# Expand the background noise to an array of the same size as the data:
	maxlenb <- max(data$lenb)
	data$bgns <- matrix(data$bgns, nrow=maxlenb, ncol=data$numb[1], byrow=TRUE)	
	
	# If the required information is present in 'data', add near range noise:
	data$bgns <- data$bgns + data$nr0a[seq_len(maxlenb),]
	
	
	##### Execution and output #####
	# Declare the lower schooling threshold, and the probability values to be returned:
	data$lsth <- zeros(dim(data$bgns), numt*lf)
	data$pr0s <- zeros(maxlenb, data$numb[1], numt*lf)
	# Define the list of indexes for the voxels for which the lower schooling threshold was set to the background noise:
	data$lst0 <- vector("list", numt*lf)
		
	# Move through the time steps:
	# Define parameters used when plotting the time bar for the generation of bottom points:
	
	infostring <- "Processing time steps:"
	cat(infostring,"\n",sep="")
	totalsteps <- length(t)
	stepfact <- nchar(infostring)/totalsteps
	oldvalue <- 0
	
	# Run through the time steps and calculate probability values:
	for(i in seq_along(t)){
		
		# Print a dot if the floor of the new value exceeds the old value in:
		thisvalue <- floor(i*stepfact)
		if(thisvalue > oldvalue){
			cat(rep(".",thisvalue-oldvalue),if(i==totalsteps) "\n", sep="")
			oldvalue <- thisvalue
		}
		
		# Estimate the expected periodic noise:
		data$pdns <- echoIBM.pdns(data,indt=i)$pdns
		
		# Define the noise of the current time step:
		#This did not work, since it reduced detectability for schools at long range: thisbgns=data$bgns * 100
		thisbgns <- data$bgns * 100
		# Add the periodic noise to the background noise:
		thisbgns[,periodic] <- thisbgns[,periodic] + data$pdns[,periodic]
		
		# Add TVG amplification to the noise:
		thisbgns <- apply.TVG(thisbgns,data,TVG.exp=TVG.exp)
				
		# Define the temporary acoustic data, used in the function, and subtract noise from the original acoustic data if required (default):
		thisvbsc <- data$vbsc[,,i]
		if(subtractNoise){
			data$vbsc[,,i] <- data$vbsc[,,i] - thisbgns
		}
		notisna <- !is.na(thisvbsc)
		p <- NAs(olddim[1:2])
		maxint <- NAs(olddim[1:2])
		
		for(f in seq_len(lf)){
			# Expand beta0 to the desired dimension:
			data$lsth[,,i+numt*(f-1)] <- beta0[i]*factor[f]
			# Wherever the background noise exceeds the lower schooling threshold, insert the background noise in data$lsth:
			data$lst0[[i+numt*(f-1)]] <- which(thisbgns>data$lsth[,,i+numt*(f-1)])
			data$lsth[,,i+numt*(f-1)][data$lst0[[i+numt*(f-1)]]] <- thisbgns[data$lst0[[i+numt*(f-1)]]]
			# Add the factor to the upper schooling threshold:
			thisbeta1 <- beta1[i]*factor[f]
			
			# The cummulative probability is 1 for data$lsth[,,i+numt*(f-1)]>=thisbeta1, which is where the background noise exceeds the upper schooling threshold (note that data$lsth[[f]] has been set to the noise due to the requirement applied with "data$lst0[[f]]" above):
			not1=notisna & data$lsth[,,i+numt*(f-1)]<thisbeta1
			is1=notisna & data$lsth[,,i+numt*(f-1)]>=thisbeta1
		
			# Run the cummulative distribution function for the normalizing constant 'maxint' and the input:
			p[not1] <- P(data$lsth[,,i+numt*(f-1)][not1],thisbeta1,thisvbsc[not1],thisbgns[not1])
			
			# Add ones where data$lsth[,,i+numt*(f-1)]>=thisbeta1:
			p[is1] <- 1
			# Add ones at the positions of high intensity noise:
			if(length(data$hini)>0){
				# Select only the high intensity noise of the current time step, and skip if the high intensity noise is not present:
				thesehini <- data$hini[data$hini[,3]==t[i],1:2,drop=FALSE]
				if(length(thesehini)>0){
					# Expand the high intensity noise along the beams:
					if(hins_add>0){
						
						newhini <- NULL
						uniquebeams <- unique(thesehini[,2])
						for(b in uniquebeams){
							thesej <- sort(thesehini[thesehini[,2]==b,1])
							diffj <- diff(thesej)
							jump <- which(diffj>1)
							if(length(jump)==0){
								lower <- 1
								upper <- length(thesej)
							}
							else{
								lower <- c(1,jump+1)
								upper <- c(jump,length(thesej))
							}
							# Expand the high intensity noise voxels along beams:
							expanded=NULL
							for(l in seq_along(lower)){
								expanded <- c(expanded, (thesej[lower[l]]-hins_add):(thesej[upper[l]]+hins_add))
							}
							# Uniquify the expanded voxels:
							expanded <- intersect(expanded, seq_len(maxlenb))
							expanded <- cbind(unique(expanded),b)
							newhini <- rbind(newhini,expanded)
						}
						thesehini <- newhini
					}
					p[arr.ind2ind(thesehini,c(maxlenb,data$numb[1]))] <- 1
				}
			}
			
			# Insert the the probability data:
			data$pr0s[,,i+numt*(f-1)] <- p
		}
	}
	
	
	# Smooth the logarithm of the probabilities (accumulated probability):
	if(length(h)>0){
		data$pr0s <- log(data$pr0s)
		# If the acoustic system is not a sonar (in which case each ping may be smoothed simultaneously using the option 'sim'), the positions data should be used several times is the length of 'factor' exceeds 1:
		if(!sim && length(factor)>1){
			sim <- length(factor)
		}
		data <- ksmooth.SG.TSD(data, ind=smind, h=h, w=3, sim=sim, drop=TRUE, var="pr0s", na.rm=na.rm)
		# Return to non-log scale:
		data$pr0s <- exp(data$pr0s)
	}
	
	# Read the acoustic data:
	if(!sim){
		temp <- read.event(event=event, cruise=cruise, t=t, var=c("voxels"), msg=FALSE, allow.old=allow.old, TOV=TOV)
		data[names(temp)] <- temp
	}
	
	# Set porbabilities outside of the subset defined by 'ind' to 1:
	ind <- ind.expand(ind,c(maxlenb, data$numb[1]))
	# These indices are mostly aimed at a specified depth range in which to run the segmentatin mehtod, and not for restictions x and y direction:
	# Update to the global coordinate system if initially read in the coordinate system of the vessel (see from line 110):
	if(cs=="v"){
		temp <- read.event(event=event, cruise=cruise, esnm=esnm, t=t, var=c("psxx", "psyx", "pszx"), cs="g", msg=FALSE, allow.old=allow.old, TOV=TOV)
		data[names(temp)] <- temp
	}
	indFromRangeSubset <- subset_TSD(data[c("psxx","psyx","pszx")], range=range, subset=subset, ind.out=TRUE)$subs[[1]]
	ind <- intersect(arr.ind2ind(ind,c(maxlenb,data$numb[1])),indFromRangeSubset)
	
	ind <- rep(ind,numt*lf) + maxlenb*data$numb[1] * rep(seq(0,numt*lf-1),each=length(ind))
	
	#plotl(c(data$pr0s[,,dim(data$pr0s)[3]]))
	data$pr0s_old <- data$pr0s
	data$pr0s <- ones(dim(data$pr0s))
	data$pr0s[ind] <- data$pr0s_old[ind]
	# Remove the temporary array:
	data$pr0s_old <- NULL
	#lines(c(data$pr0s[,,dim(data$pr0s)[3]]), col=2)
	
	# Apply the threshold:
	if(length(alpha)>0){
		data$sgsc <- which(data$pr0s<alpha)
	}
	
	# Return the output:
	data
	##################################################
	##################################################
}

echoIBM.vbsc2p.event_old <- function(event=1, t=1, cruise=2009116, bgns=NULL, beta0=NULL, beta1=NULL, factor=NULL, ind=list(-(1:30),NULL), nsind=0.75, smind=list(-(1:300)), range=list(), subset=NULL, h=NULL, alpha=NULL, dir.data=NULL, hins_add=10, phase=TRUE, TVG.exp=2, esnm="MS70", subtractNoise=TRUE, adds=list(), sim=TRUE, na.rm=1, allow.old=FALSE, TOV=0){

	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-05-23 - Clean version.
	# Update: 2012-06-03 - Added a solution to the problem of missing values due to the limited range over which the exponential integral calculated by expint_E1() returns values with accepted error. This range is [-701.8114,701.8114]. For values of the argument in the second exponential integral expression, -thisx*(thisbeta0-thisbetaN)/thisbeta0/thisbetaN < -700, 'betaN', the noise level, is adjusted so that the value at -700 is used. This is done as an approsimation, after the observation that the effect of the background noise is far weaker than the effect of the signal x. This approximation will allways give a higher probability that the signal<beta, which represents a conservative estimate. Thus the school is less likely to be detected in the voxels for which the approximation is applied. However, this happens for low noise levels and high acoustic values, which suggests that the school is detected anyhow. See the script "echoIBM.vbsc2p.event_approximation.R" for plot of the approximation.
	# Update: 2012-06-03 - Changed the calculation to the new expression of the probability used in paper 2.
	# Update: 2012-07-02 - Added high intensity noise.
	# Update: 2012-07-12 - Added the option of expanding the high intensity noise along beams by 'hins_add' voxels, to include high values that are not high enough to be detected as high intensity noise values.
	# Update: 2012-05-23 - Added the parameter 'hins_add' for discarding voxels on both sides along beams of the inentified high intensity noise.
	# Update: 2012-05-23 - Added support for multiple time steps.
	# Update: 2012-11-12 - Restructured to estimate the phase at each ping.
	# Update: 2012-11-14 - Substituted the option 'pns3new' by 'phase' with default TRUE, indicating estimation of the phase at each ping.
	# Update: 2012-12-18 - Added 'insert.NA' and 'set2noise'.
	# Update: 2013-01-03 - Removed 'insert.NA' and 'set2noise'.
	# Update: 2013-05-09 - Added support for subsetting by 'ind', 'range', and 'subset'.
	# Update: 2013-08-12 - Added the option 'factor', useful for calculating multiple segmentation masks at once.
	# Update: 2013-09-12 - Added the option 'sim'.
	# Last: 2015-06-03 - Added center of mass in the output.


	##################################################
	##################################################
	##### Preparation #####
	# A small funciton used to get the minimum of each row in a two column matrix (much faster than apply(x,1,min)):
	rowMin2 <- function(x){
		a = x[,2] > x[,1]
		x[a,2] = x[a,1]
		x[,2]
		}

	# Define the function for calculating the cummulative distribution. This function takes up to one second for one ping of MS70 data, mostly due to the expint_E1() function:
	P <- function(thisbeta0, thisbeta1, thisx, thisbetaN){
		# Define the maximum value at which the function expint_E1() returns non-missing values:
		maxarg = 700
		# Identify the values giving missing values in expint_E1() for 'beta0' and 'beta1':
		outsidemaxarg_0 = thisx*(thisbeta0-thisbetaN)/thisbeta0/thisbetaN > maxarg
		outsidemaxarg_1 = thisx*(thisbeta1-thisbetaN)/thisbeta1/thisbetaN > maxarg

		# Adjust the background noise to reach non-missing values
		thisbetaN[outsidemaxarg_0] = 1/(1/thisbeta0[outsidemaxarg_0] + maxarg/thisx[outsidemaxarg_0])
		thisbetaN[outsidemaxarg_1] = 1/(1/thisbeta1 + maxarg/thisx[outsidemaxarg_1])
		# Define the delta values used in the formula:
		deltabeta0 = 1/thisbeta0 - 1/thisbetaN
		deltabeta1 = 1/thisbeta1 - 1/thisbetaN

		# If absolute values still are outside of 'maxarg', or the argument to the first expint_E1() is outside of the definition area, insert .Machine$double.eps at those positions:
		outsidemaxarg_0 = thisx*(1/thisbeta0-1/thisbetaN) > maxarg | thisx*(1/thisbeta1-1/thisbetaN) > maxarg | thisx/thisbeta0 > maxarg | thisx/thisbeta1 > maxarg

		H0 = expint_E1(thisx/thisbeta0) - exp(-thisx/thisbetaN) * (expint_E1(thisx*deltabeta0) + log(abs(thisbeta0*deltabeta0)))
		H1 = expint_E1(thisx/thisbeta1) - exp(-thisx/thisbetaN) * (expint_E1(thisx*deltabeta1) + log(abs(thisbeta1*deltabeta1)))

		out = H0/H1
		out[outsidemaxarg_0] = .Machine$double.eps
		out
		}


	lf = length(factor)
	numt = length(t)

	# 'beta0' or 'beta1' should be of the same length as the number of time steps:
	if(length(beta0)!=numt){
		beta0 = rep(beta0,length.out=numt)
		}
	if(length(beta1)!=numt){
		beta1 = rep(beta1,length.out=numt)
		}


	# Update 'esnm':
	esnm = read.event(event=event, cruise=cruise, esnm=esnm, var="esnm", msg=FALSE)$esnm

	# 'sim' can only be TRUE for sonars:
	if(!is.sonar(esnm)){
		sim=FALSE
		}

	# For simultaneous smoothing for all time steps using the voxels positions in the coordinate system of the vessel, read the voxels positions only once:
	#if(sim){
	#	dynamicvar=c("vbsc","time")
	#	staticvar=c("psxx","psyx","pszx","volx","beams")
	#	cs="v"
	#	}
	#else{
	#	dynamicvar=c("vbsc","time","esnm","psxx","psyx","pszx")
	#	staticvar=c("volx","beams")
	#	cs="g"
	#	}
	if(sim){
		dynamicvar=c("vbsc", "time", "beams")
		staticvar=c("psxx", "psyx", "pszx", "volx")
		cs="v"
		}
	else{
		dynamicvar=c("vbsc", "time", "esnm", "psxx", "psyx", "pszx", "beams")
		staticvar=c("volx")
		cs="g"
		}
	noisevar=c("hini","hins")

	# Read the acoustic data:
	data = read.event(event=event, cruise=cruise, esnm=esnm, t=t, var=dynamicvar, msg=FALSE, cs=cs, allow.old=allow.old, TOV=TOV)
	# Add high intensity noise:
	temp = read.event(event=event, cruise=cruise, t=t, var=noisevar, msg=FALSE, allow.old=allow.old)
	data[names(temp)] = temp
	# Read the voxel position data in the coordinate system of the vessel, which is static for all time steps, allowing for the smoothing to be done over several time steps at the time. Also add the sonar configuration:
	temp = read.event(event=event, cruise=cruise, t=1, var=staticvar, msg=FALSE, cs=cs, allow.old=allow.old)
	data[names(temp)] = temp
	# Store the original dimension of the acoustic data:
	if(length(dim(data$vbsc))==2){
		dim(data$vbsc) = c(dim(data$vbsc),numt)
		}
	olddim = dim(data$vbsc)
	# Read background and periodic noise data, added to the optional list 'bgns':
	noisefiles = noise.path(cruise=cruise,event=event,esnm=esnm,dir.data=dir.data)
	data = c(adds, data, bgns, read.TSDs(noisefiles,var=c("bgns","badb","pns1","pns2","pns3","pn3M","harm","acfq","nr0a")))
	periodic = which(data$badb==1)
	rm(bgns)
	gc()

	# If zeros occur in data$bgns, add a small value:
	if(min(data$bgns)==0){
		data$bgns[data$bgns==0] = .Machine$double.eps
		}

	# If the required information is present in 'data', estimate the phase of the periodic noise at each time step, from a portion of the sonar volume dominated by noise, or alternatively, not occupied by schools:
	if(phase){
		data$pn3M = get.pdns_phase.event(event=event, cruise=cruise, esnm=esnm, t=t, noise=data, nsind=nsind)$pn3M
		}

	# Expand the background noise to an array of the same size as the data:
	maxlenb = max(data$lenb)
	data$bgns = matrix(data$bgns, nrow=maxlenb, ncol=data$numb[1], byrow=TRUE)	

	# If the required information is present in 'data', add near range noise:
	data$bgns = data$bgns + data$nr0a[seq_len(maxlenb),]


	##### Execution and output #####
	# Declare the lower schooling threshold, and the probability values to be returned:
	data$lsth = zeros(dim(data$bgns), numt*lf)
	data$pr0s = zeros(maxlenb, data$numb[1], numt*lf)
	# Define the list of indexes for the voxels for which the lower schooling threshold was set to the background noise:
	data$lst0 = vector("list", numt*lf)

	# Move through the time steps:
	# Define parameters used when plotting the time bar for the generation of bottom points:

	infostring = "Processing time steps:"
	cat(infostring,"\n",sep="")
	totalsteps = length(t)
	stepfact = nchar(infostring)/totalsteps
	oldvalue = 0

	# Run through the time steps and calculate probability values:
	for(i in seq_along(t)){

		# Print a dot if the floor of the new value exceeds the old value in:
		thisvalue = floor(i*stepfact)
		if(thisvalue > oldvalue){
			cat(rep(".",thisvalue-oldvalue),if(i==totalsteps) "\n", sep="")
			oldvalue = thisvalue
			}

		# Estimate the expected periodic noise:
		data$pdns = echoIBM.pdns(data,indt=i)$pdns

		# Define the noise of the current time step:
		#This did not work, since it reduced detectability for schools at long range: thisbgns=data$bgns * 100
		thisbgns = data$bgns * 100
		# Add the periodic noise to the background noise:
		thisbgns[,periodic] = thisbgns[,periodic] + data$pdns[,periodic]

		# Add TVG amplification to the noise:
		thisbgns = apply.TVG(thisbgns,data,TVG.exp=TVG.exp)
		
		# Define the temporary acoustic data, used in the function, and subtract noise from the original acoustic data if required (default):
		thisvbsc = data$vbsc[,,i]
		if(subtractNoise){
			data$vbsc[,,i] = data$vbsc[,,i] - thisbgns
			}
		notisna = !is.na(thisvbsc)
		p = NAs(olddim[1:2])
		maxint = NAs(olddim[1:2])

		for(f in seq_len(lf)){
			# Expand beta0 to the desired dimension:
			data$lsth[,,i+numt*(f-1)] = beta0[i]*factor[f]
			# Wherever the background noise exceeds the lower schooling threshold, insert the background noise in data$lsth:
			data$lst0[[i+numt*(f-1)]] = which(thisbgns>data$lsth[,,i+numt*(f-1)])
			data$lsth[,,i+numt*(f-1)][data$lst0[[i+numt*(f-1)]]] = thisbgns[data$lst0[[i+numt*(f-1)]]]
			# Add the factor to the upper schooling threshold:
			thisbeta1 = beta1[i]*factor[f]
	
			# The cummulative probability is 1 for data$lsth[,,i+numt*(f-1)]>=thisbeta1, which is where the background noise exceeds the upper schooling threshold (note that data$lsth[[f]] has been set to the noise due to the requirement applied with "data$lst0[[f]]" above):
			not1=notisna & data$lsth[,,i+numt*(f-1)]<thisbeta1
			is1=notisna & data$lsth[,,i+numt*(f-1)]>=thisbeta1

			# Run the cummulative distribution function for the normalizing constant 'maxint' and the input:
			p[not1] = P(data$lsth[,,i+numt*(f-1)][not1],thisbeta1,thisvbsc[not1],thisbgns[not1])
	
			# If there are pairs of 'data$lsth[,,i+numt*(f-1)]' and 'thisbgns' that are identical, use a simple average of the closest value to either side of 'thisbgns' (Cauchy principal value):
			equal = which(data$lsth[,,i+numt*(f-1)][not1]==thisbgns[not1])
			if(length(equal)){
				# Get the small value to shift the estimation of the probabilities with to both sides in the case that the lower schooling threshold has been set to the background noise:
				add = rowMin2(cbind(thisbgns[not1][equal]/2,sqrt(.Machine$double.eps)))
				# Apply the Cauchy principal value:
				p[not1][equal] = rowMeans( cbind(P(thisbgns[not1][equal]-add,thisbeta1,thisvbsc[not1][equal],thisbgns[not1][equal]), P(thisbgns[not1][equal]+add,thisbeta1,thisvbsc[not1][equal],thisbgns[not1][equal])), na.rm=TRUE)
				#p[not1][equal]=rowMeans( cbind(P(thisbgns[not1][equal]-sqrt(.Machine$double.eps),thisbeta1,thisvbsc[not1][equal],thisbgns[not1][equal]) + P(thisbgns[not1][equal]+sqrt(.Machine$double.eps),thisbeta1,thisvbsc[not1][equal],thisbgns[not1][equal])), na.rm=TRUE)
				}
	
			# Add ones where data$lsth[,,i+numt*(f-1)]>=thisbeta1:
			p[is1] = 1
			# Add ones at the positions of high intensity noise:
			if(length(data$hini)>0){
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
							expanded=NULL
							for(l in seq_along(lower)){
								expanded = c(expanded, (thesej[lower[l]]-hins_add):(thesej[upper[l]]+hins_add))
								}
							# Uniquify the expanded voxels:
							expanded = intersect(expanded, seq_len(maxlenb))
							expanded = cbind(unique(expanded),b)
							newhini = rbind(newhini,expanded)
							}
						thesehini = newhini
						}
					p[arr.ind2ind(thesehini,c(maxlenb,data$numb[1]))] = 1
					}
				}
	
			# Insert the the probability data:
			data$pr0s[,,i+numt*(f-1)] = p
			}
		}


	# Smooth the logarithm of the probabilities (accumulated probability):
	if(length(h)>0){
		data$pr0s = log(data$pr0s)
		# If the acoustic system is not a sonar (in which case each ping may be smoothed simultaneously using the option 'sim'), the positions data should be used several times is the length of 'factor' exceeds 1:
		if(!sim && length(factor)>1){
			sim = length(factor)
			}
		data = ksmooth.SG.TSD(data, ind=smind, h=h, w=3, sim=sim, drop=TRUE, var="pr0s", na.rm=na.rm)
		# Return to non-log scale:
		data$pr0s = exp(data$pr0s)
		}

	# Read the acoustic data:
	if(!sim){
		temp = read.event(event=event, cruise=cruise, t=t, var=c("voxels"), msg=FALSE, allow.old=allow.old, TOV=TOV)
		data[names(temp)] = temp
		}

	# Set porbabilities outside of the subset defined by 'ind' to 1:
	ind = ind.expand(ind,c(maxlenb, data$numb[1]))
	# These indices are mostly aimed at a specified depth range in which to run the segmentatin mehtod, and not for restictions x and y direction:
	# Update to the global coordinate system if initially read in the coordinate system of the vessel (see from line 110):
	if(cs=="v"){
		temp = read.event(event=event, cruise=cruise, esnm=esnm, t=t, var=c("psxx", "psyx", "pszx"), cs="g", msg=FALSE, allow.old=allow.old, TOV=TOV)
		data[names(temp)] = temp
		}
	indFromRangeSubset = subset_TSD(data[c("psxx","psyx","pszx")], range=range, subset=subset, ind.out=TRUE)$subs[[1]]
	ind = intersect(arr.ind2ind(ind,c(maxlenb,data$numb[1])),indFromRangeSubset)

	ind = rep(ind,numt*lf) + maxlenb*data$numb[1] * rep(seq(0,numt*lf-1),each=length(ind))

	#plotl(c(data$pr0s[,,dim(data$pr0s)[3]]))
	data$pr0s_old = data$pr0s
	data$pr0s = ones(dim(data$pr0s))
	data$pr0s[ind] = data$pr0s_old[ind]
	# Remove the temporary array:
	data$pr0s_old = NULL
	#lines(c(data$pr0s[,,dim(data$pr0s)[3]]), col=2)

	# Apply the threshold:
	if(length(alpha)>0){
		data$sgsc = which(data$pr0s<alpha)
		}

	# Return the output:
	data
}
