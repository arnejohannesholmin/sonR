#*********************************************
#*********************************************
#' Locates the point at which the kernel desity estimate of the Sv of the subset is equal to half its maximum, and uses Paper III of the PhD of Holmin to estimate the mean Sv from this. The output variables all start with X, denoting the SX90 segmentation method, but are applicable to all segmentation methods.
#'
#' @param data			A list containing the acoustic data eihter given as 'vbss' for a segment of the data or as 'vbsc'.
#' @param plot.hist		Logical: If TRUE plot the histogram of the Sv.
#' @param allow.vbsc	?
#' @param list.out		Use a list as output.
#' @param minlen		The minimum length of the data, at and below which NA is returned.
#' @param type			A single character denoting the type of segmentation to label the output with.
#' @param enlarged		Logical: If TRUE add "E" to the variable names and remove the last character.
#' @param ...			Passed on to \code{\link{density}} and \code{\link{hist}}.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom stats median approx optimize quantile uniroot
#' @importFrom graphics lines
#'
#' @export
#' @rdname meanSv.TSD
#'
meanSv.TSD <- function(data, plot.hist=FALSE, allow.vbsc=TRUE, list.out=FALSE, minlen=1, type="H", enlarged=FALSE, ...){
	
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
	##### Preparation #####
	# The scale of the negative Gumbel distribution fitting the Sv:
	sigma=10/log(10)
	
	# Calculate the kernel density estimate:
	if(length(data$vbss)>1){
		suppressWarnings(x <- 10*log10(data$vbss))
		if(length(x)>1){
			suppressWarnings(kd<-density(x, na.rm=TRUE, ...))
			}
		else{
			kd=x
			}
		}
	else if(length(data$sgsc)>1 && length(data$vbsc)>1){
		suppressWarnings(x <- 10*log10(data$vbsc[data$sgsc]))
		if(length(x)>1){
			suppressWarnings(kd<-density(x, na.rm=TRUE, ...))
			}
		else{
			kd=x
			}
		}
	else if(allow.vbsc && length(data$vbsc)>1){
		suppressWarnings(x <- 10*log10(data$vbsc))
		if(length(x)>1){
			suppressWarnings(kd<-density(x, na.rm=TRUE, ...))
			}
		else{
			kd=x
			}
		}
	else if(allow.vbsc && length(data$mvbs)>1){
		suppressWarnings(x <- data$mvbs)
		if(length(x)>1){
			suppressWarnings(kd<-density(x, na.rm=TRUE, ...))
			}
		else{
			kd=x
			}
		}
	#else if(length(data)>0){
	#	data=list(vbss=10*log10(data))
	#	}
	else{
		warning("No voxels segmented")
		return(list())
		}
	if(plot.hist){
		ll=list(...)
		thisl=list(x=x,plot=FALSE)
		otherl=ll[setdiff(names(ll),names(thisl))]
		hh=do.call("hist",c(thisl,otherl))
		thisl=list(x=x,probability=TRUE,ylim=range(hh$density,kd$y,na.rm=TRUE),axes=TRUE,xlab="x",ylab="Probability",main="Histogram")
		otherl=ll[setdiff(names(ll),names(thisl))]
		hh=do.call("hist",c(thisl,otherl))
		lines(kd)
		}
	
	##### Execution and output #####
	# Locate the maximum value and mirror the values to the left of the peak across this value:
	atmaxkd=which.max(kd$y)
	maxkd=max(kd$y)
	# Use this to discard the lowes half point:
	kd$w = kd$y
	kd$w[seq_len(atmaxkd)]=2*maxkd-kd$y[seq_len(atmaxkd)]
	# Create a function for locating the half value:
	f=function(x){
		approx(kd$x, kd$w, x)$y - 0.5 * maxkd
		}
	# Locate the value of Sv at the half value:
	Svhalf=uniroot(f,c(kd$x[atmaxkd],max(kd$x)))$root
	
	# Function that calculates f(Svhalf)-exp(-1)/2:
	halfXpSvGumbelSv=function(beta,Svhalf){
		mu=sigma*log(beta)
		frac=(Svhalf-mu)/sigma
		abs( log( log(2*exp(1)) + frac ) - frac )
		}
	svhalf=10^(Svhalf/10)
	k=1e3
	o=optimize(halfXpSvGumbelSv,interval=c(svhalf/k,svhalf),Svhalf=Svhalf)
	
	HhSv = 10*log10(o$m)
	
	# Fast approximation, where the scale nu = 10/log(10), a quantile of the kernel density estimate is taken at 93 percent (giving a point close to the half point for a negative Gumbel distributed variable), and the Gumbel CDF is used to calculate the mean mu, by GumbelCDF(y at the 93 percentile) = exp( -exp(-(y-mu)/nu) ) = 0.07 => y-mu = 10/log(10) * log(-log(0.07)) = 4.247608 approx 4.25. Thus, the mean Sv can be estimated by the 93-percentile of the Sv minus 4.25 dB:
	HGSv = quantile(x, 0.93, names=FALSE, na.rm=TRUE) - 4.25
	
	# Also get the average and median Sv:
	HaSv = 10*log10(mean(10^(x/10), na.rm=TRUE))
	HdSv = 10*log10(median(10^(x/10), na.rm=TRUE))
	
	# Get the quantile of this point:
	#Qkern = sum(kd$y[kd$x<=Svhalf]) / sum(kd$y)
	#Q = mean(x<=Svhalf)
	
	# Get the mean Sv from the peak, which is at nu * log(beta):
	HpSv = 10*log10(exp(kd$x[atmaxkd] / sigma))
		
	# Also calculate the convensional mean Sv:
	if(length(data$vbsc)>0 && length(data$volx)>0){
		HmSv=10*log10(sum(data$vbsc[data$sgsc] * data$volx[data$sgsc], na.rm=TRUE) / sum(data$volx[data$sgsc], na.rm=TRUE))
		#HmSv=NULL
		}
	else if(length(data$vbss)>0 && length(data$vols)>0){
		HmSv=10*log10(sum(data$vbss*data$vols,na.rm=TRUE)/sum(data$vols,na.rm=TRUE))
		#HmSv=NULL
		}
	else{
		HmSv=NA
		}
		
	# If fewer points than the specified minimum number of points is present, return NA for the variables based on the kernel density estimate:
	if(length(x)<minlen){
		HpSv = NA
		HhSv = NA
		}
	
	# Convert to linear values:	
	Hasv = 10^(HaSv/10) # Average
	Hdsv = 10^(HdSv/10) # Median
	Hmsv = 10^(HmSv/10) # Mean
	Hpsv = 10^(HpSv/10) # Peak
	Hhsv = 10^(HhSv/10) # Half peak
	HGsv = 10^(HGSv/10) # Gumbel quantile
	

	# Return:
	if(list.out){
		out = list(Hasv=Hasv, HaSv=HaSv, Hdsv=Hdsv, HdSv=HdSv, Hmsv=Hmsv, HmSv=HmSv, Hpsv=Hpsv, HpSv=HpSv, Hhsv=Hhsv, HhSv=HhSv, HGsv=HGsv, HGSv=HGSv)
		}
	else{
		out = c(Hasv=Hasv, HaSv=HaSv, Hdsv=Hdsv, HdSv=HdSv, Hmsv=Hmsv, HmSv=HmSv, Hpsv=Hpsv, HpSv=HpSv, Hhsv=Hhsv, HhSv=HhSv, HGsv=HGsv, HGSv=HGSv)
		}
	# Add names:
	namesend = c("asv", "aSv", "dsv", "dSv", "msv", "mSv", "psv", "pSv", "hsv", "hSv", "Gsv", "GSv")
	namesout = paste0(type[1], if(enlarged) "E", if(enlarged) substr(namesend,1,2) else namesend)
	names(out) <- namesout
	out
	##################################################
	##################################################
	}
