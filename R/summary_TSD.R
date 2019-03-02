#*********************************************
#*********************************************
#' Calculates variaous measures of mean sv, and total backscatter and volume.
#'
#' @param data		A list containing volume and the acoustic data eihter given as 'vbss' for a segment of the data or as 'vbsc'.
#' @param plot.hist	Logical: If TRUE plot the histogram when calculating the robust mean sv.
#' @param minlen	The minimum number of acoustic samples.
#' @param ...		Passed on to density() when calculating the robust mean sv.
#'
#' @importFrom data.table rbindlist
#' @importFrom stats median approx density optimize quantile sd uniroot
#' @importFrom graphics lines
#'
#' @export
#'
summary_TSD <- function(data, plot.hist=FALSE, minlen=1, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-01-02 - Clean version.
	# Last: 2015-07-01 - Added all possible estimation methods, and assigned TSD names to the output:
	########### DESCRIPTION: ###########
	# Locates the point at which the kernel desity estimate of the Sv of the subset is equal to half its maximum, and uses Paper III of the PhD of Holmin to estimate the mean Sv from this. The output variables all start with X, denoting the SX90 segmentation method, but are applicable to all segmentation methods.
	########## DEPENDENCIES: ###########
	# 
	############ VARIABLES: ############
	# ---data--- is a list containing the acoustic data eihter given as 'vbss' for a segment of the data or as 'vbsc'.
		
	
	##################################################
	##################################################
	# Funciton for extracting mean sv, total sv and total volume of one ping:
	getSummaryOnePing <- function(data, t, sigma, plot.hist=FALSE, minlen=1, ...){
		##### Preparation #####
		# Extract the requested ping:
		data <- extractTimeStep(data, t=t)
		
		# Extract the volume of the voxels:
		volx <- data$volx
		harx <- data$harx
		
		# Extract the volume backscattering data:
		if(length(data$vbss)>1){
			vbsc <- data$vbss
			volx <- NA
			harx <- NA
		}
		else if(length(data$sgsc)>1 && length(data$vbsc)>1){
			vbsc <- data$vbsc[data$sgsc]
			volx <- volx[data$sgsc]
			harx <- harx[data$sgsc]
		}
		else if(length(data$vbsc)>1){
			vbsc <- data$vbsc
		}
		else if(length(data$mvbs)>1){
			vbsc<- 10^(data$mvbs/10)
		}
		else{
			warning("Acoustic data missing")
			return(list())
		}
		
		# Set volumes at missing volume backscattering data to NA:
		volx[is.na(data$vbsc)] <- NA
		harx[is.na(data$vbsc)] <- NA
		
		# Use decebel data for the kernel stuff:
		mvbs <- 10 * log10(vbsc)
		
		# Calculate the kernel density estimate:
		if(length(mvbs)>1){
			suppressWarnings(kd <- density(mvbs, na.rm=TRUE, ...))
			}
		else{
			kd <- mvbs
			}
		if(plot.hist){
			ll <- list(...)
			thisl <- list(x=kd, plot=FALSE)
			otherl <- ll[setdiff(names(ll), names(thisl))]
			hh <- do.call("hist",c(thisl,otherl))
			thisl <- list(x=kd, probability=TRUE, ylim=range(hh$density,kd$y,na.rm=TRUE), axes=TRUE, xlab="x", ylab="Probability", main="Histogram")
			otherl <- ll[setdiff(names(ll), names(thisl))]
			hh <- do.call("hist", c(thisl,otherl))
			lines(kd)
			}
	
	
		##### Execution and output #####
		# Locate the maximum value and mirror the values to the left of the peak across this value:
		atmaxkd <- which.max(kd$y)
		maxkd <- max(kd$y)
		# Use this to discard the lowes half point:
		kd$w <- kd$y
		kd$w[seq_len(atmaxkd)] <- 2 * maxkd - kd$y[seq_len(atmaxkd)]
		# Create a function for locating the half value:
		f <- function(x){
			approx(kd$x, kd$w, x)$y - 0.5 * maxkd
			}
		# Locate the value of Sv at the half value:
		Svhalf <- uniroot(f, c(kd$x[atmaxkd], max(kd$x)))$root
	
		# Function that calculates f(Svhalf)-exp(-1)/2:
		halfXpSvGumbelSv <- function(beta, Svhalf){
			mu <- sigma * log(beta)
			frac <- (Svhalf - mu) / sigma
			abs( log( log(2 * exp(1)) + frac ) - frac )
			}
		svhalf <- 10^(Svhalf/10)
		k <- 1e3
		o <- optimize(halfXpSvGumbelSv, interval=c(svhalf/k,svhalf), Svhalf=Svhalf)
	
		rbsv <- o$m
	
		# Fast approximation, where the scale nu = 10/log(10), a quantile of the kernel density estimate is taken at 93 % (giving a point close to the half point for a negative Gumbel distributed variable), and the Gumbel CDF is used to calculate the mean mu, by GumbelCDF(y at the 93 percentile) = exp( -exp(-(y-mu)/nu) ) = 0.07 => y-mu = 10/log(10) * log(-log(0.07)) = 4.247608 approx 4.25. Thus, the mean Sv can be estimated by the 93-percentile of the Sv minus 4.25 dB:
		Gmsv <- quantile(mvbs, 0.93, names=FALSE, na.rm=TRUE) - 4.25
	
		# Also get the average and median Sv:
		avsv = mean(vbsc, na.rm=TRUE)
		mdsv = median(vbsc, na.rm=TRUE)
		mxsv = max(vbsc, na.rm=TRUE)
	
		# Get the quantile of this point:
		#Qkern = sum(kd$y[kd$x<=Svhalf]) / sum(kd$y)
		#Q = mean(x<=Svhalf)
	
		# Get the mean Sv from the peak, which is at nu * log(beta):
		pksv <- exp(kd$x[atmaxkd] / sigma)
		
		# Also calculate the convensional mean Sv:
		tbsc <- sum(vbsc * volx, na.rm=TRUE)
		tbsc_harx <- sum(vbsc * harx, na.rm=TRUE)
		tvol <- sum(volx, na.rm=TRUE)
		thar <- sum(harx, na.rm=TRUE)
		mesv <- tbsc / tvol
		masv <- tbsc_harx / thar
		
		sdsv <- sd(vbsc, na.rm=TRUE)
		sdSv <- sd(mvbs, na.rm=TRUE)
		
		# If fewer points than the specified minimum number of points is present, return NA for the variables based on the kernel density estimate:
		if(length(vbsc) < minlen){
			pksv <- NA
			hpsv <- NA
			}
	
	
		# Convert to lobarithmic values:
		avSv = 10 * log10(avsv) # Average Sv
		mdSv = 10 * log10(mdsv) # Median Sv
		meSv = 10 * log10(mesv) # Mean Sv, by sum(sv * vol) / sum(vol)
		maSv = 10 * log10(masv) # Mean Sv given the horizontal area
		rbSv = 10 * log10(rbsv) # Parametric robust Gumbel mean Sv
		pkSv = 10 * log10(pksv) # Peak Gumbel mean Sv
		GmSv = 10 * log10(Gmsv) # Gumbel approximation of mean Sv
		mxSv = 10 * log10(mxsv) # Maximum Sv
		tgTS = 10 * log10(tbsc) # Total backascatter
		# sdsv # Standard deviation of the Sv
		# sdSv # Standard deviation of the Sv
		# tvol # Total volume
		# thar # Total horizontal area
	

		# Return:
		out = data.frame(
			avSv=avSv, avsv=avsv,	
			mdSv=mdSv, mdsv=mdsv,	
			meSv=meSv, mesv=mesv,	
			maSv=maSv, masv=masv,	
			rbSv=rbSv, rbsv=rbsv,	
			pkSv=pkSv, pksv=pksv,	
			GmSv=GmSv, Gmsv=Gmsv,	
			mxSv=mxSv, mxsv=mxsv,	
			sdSv=sdSv, sdsv=sdsv,	
			tgTS=tgTS, tbsc=tbsc,	
			tvol=tvol,
			thar=thar)
						 
		return(out)
	}
	
	# The scale of the negative Gumbel distribution fitting the Sv:
	sigma=10/log(10)
	
	# Get the number of time steps:
	if(length(data$numt)){
		numt <- data$numt
	}
	else if(length(data$utim)){
		numt <- length(data$utim)
	}
	else{
		warning("'data' needs either 'numt' or 'utim'. Number of time steps set to 1")
		numt <- 1
	}	


	# Run through the time steps:
	s <- lapply(seq_len(numt), function(i) getSummaryOnePing(data, t=i, sigma=sigma, plot.hist=plot.hist, minlen=minlen, ...))
	data.table::rbindlist(s)
	##################################################
	##################################################
	}
