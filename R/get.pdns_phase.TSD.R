#*********************************************
#*********************************************
#' Estimates the phase of the periodic noise of the MS70 sonar. If the variables 'badb', 'pns1', 'pns2', 'harm' and 'bgns' are not present in the input list 'data', an attempt is made to read them from the directory holding the MS70 noise estimate for the given event.
#'
#' @param data  is a list of the following variables ( can be obtained by data=c(read.event(var=c("beams","vbsc","")),read.TSDs(noise.path(),var=c("pns1","pns2","harm","badb","bgns","acfq"))) ):
#' @param nsind  is a vector of indexes along the beams, as input to ind.expand(), used to select the subset over which the estimation of the phase of the periodic noise is done. If given as a single numeric, the outermost 'nsind' voxels are used in each beam.
#' @param cruise,event,esnm  Used in \code{\link{noise.path}} to get the path to the noise data holding the background noise, if missing in the data.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param pdns_scale  is used to scale the noise in order to allow the optimization to work.
#' @param TVG  is FALSE if TVG compensation is to be removed from the data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw apply.TVG
#' @importFrom TSD ind.expand NAs read.TSDs
#' @importFrom utils tail
#' @importFrom stats optimize
#'
#' @export
#' @rdname get.pdns_phase.TSD
#'
get.pdns_phase.TSD<-function(data=list(), nsind=0.75, cruise=2009116, event=NULL, esnm="MS70", dir.data=NULL, pdns_scale=1e-14, TVG=TRUE){
	
	############### LOG: ###############
	# Start: 2010-11-23 - Clean version, adopted from get.pdns_phase.MS70().
	# Update: 2013-07-19 - Changed 'ind1' to 'nsind' and applied ind.expand().
	# Last: 2013-07-19 - Changed names from get.pdns_phase.MS70.TSD() to get.pdns_phase.TSD().
	
	
	browser()
	##### Preparation #####
	# Function for calculating the mean og angles:
	angleMean=function(ang,w=1){
		atan2(sum(w*sin(ang),na.rm=TRUE),sum(w*cos(ang),na.rm=TRUE)) %% (2*pi)
		}
	
	# Read background and periodic noise data if missing in 'data':
	if(any(length(data$badb)==0, length(data$pns1)==0, length(data$pns2)==0, length(data$harm)==0, length(data$bgns)==0)){
		noisefiles=noise.path(cruise=cruise,event=event,dir.data=dir.data,esnm=esnm)
		noise=read.TSDs(noisefiles,var=c("acfq","badb","pns1","pns2","harm","bgns"),dimension=TRUE)
		data[names(noise)]=noise
		}
	# Check whether the required variables are present:
	requiredVariables=c("vbsc","acfq","badb","pns1","pns2","harm","bgns")
	if(!all(requiredVariables %in% names(data))){
		warning(paste("At least one of the following variables are missing in the input data: Volume backscattering coefficient 'vbsc', and periodic noise variables 'acfq', 'badb', 'pns1', 'pns2', 'harm', 'bgns'. Variables present are: ",paste(sort(intersect(names(data),requiredVariables)),collapse=", ") ,sep=""))
		return(list())
		}
	
	# Add one dimension to the acoustic data for convenience:
	if(length(dim(data$vbsc))==2){
		dim(data$vbsc)=c(dim(data$vbsc),1)
		}
	ntimesteps=dim(data$vbsc)[3]
	numf = length(unique(c(data$freq)))
	dimb = c(data$numb[1]/numf, numf)
	
	
	# Define the factor in the sine wave (sin(ax) gives period dt = 2*pi/a => a = 2*pi/dt. dt = time between peaks 1/acfq divided by pulselength 'sint', and divided by 2 to fit to the observed frequencies => a = 2*pi * acfq*sint):
	twopif = 2*pi * data$acfq[1] * data$sint[1]
	
	# Define the function that calculates the periodic noise added the background noise:
	pulsenoise_phase = function(phase, magnitude, kurtosis, t, n, twopif, bgns){
		bgns + magnitude*10^(kurtosis*sin(n*twopif*t + phase))
		}
	
	# Define the function used in the opimazition (bgns + pns1 * exp(pns2 * sin(n*a*j - pns3))):
	pulsenoiseSum = function(phase, x, magnitude, kurtosis, t, n, twopif, bgns){
		#sum(abs(pulsenoise(par,t,n,twopif,bgns)-x)) No longer in use
		sum((pulsenoise_phase(phase,magnitude,kurtosis,t,n,twopif,bgns)-x)^2)
		}
	

	# Define parameter interval for the optimization:
	interval = c(-2*pi, 4*pi)

	# Define the valid voxels:
	nsind = ind.expand(nsind, max(data$lenb), drop=TRUE)
	
	
	##### Execution and output #####
	pn3M = NAs(data$numb[1], ntimesteps)
	# Scale the magnitude parameter of the periodic noise according to the scaling of the acoustic data, to assure successful optimization:
	data$pns1 = data$pns1 / pdns_scale
	
	### 3. Estimate the parameters of the periodic noise: ###
	if(sum(data$badb)>0){
		periodic = which(data$badb==1)
		p = 1
		for(i in seq_len(ntimesteps)){
			# Select the acoustic data of the time step, and ignore near-range and apply the subset specified by 'nsind':
			if(TVG){
				thisvbsc = apply.TVG(data$vbsc[,,i], data, rm=TRUE)
				}
			#thisvbsc=data$vbsc[,,i]
			thisvalidvoxels = intersect(which(!apply(thisvbsc,1,function(x) any(is.na(x)))), nsind)
			if(length(pdns_scale)==0 && i==data$indt[1]){
				pdns_scale = mean(thisvbsc[thisvalidvoxels,])
				}
			thisvbsc = thisvbsc[thisvalidvoxels,,drop=FALSE] / pdns_scale
			
			# Optimize each of the beams identified to be affected by the periodic noise:
			for(j in seq_along(periodic)){
				o = optimize(f=pulsenoiseSum, interval=interval, x=thisvbsc[,periodic[j]], magnitude=data$pns1[periodic[j]], kurtosis=data$pns2[periodic[j]], t=thisvalidvoxels, n=data$harm[periodic[j]], twopif=twopif, bgns=data$bgns[periodic[j]])
				pn3M[periodic[j],p] = o$minimum %% (2*pi)
				}
			p = p+1
			}
		}
	
	# Equal the phase estimates for each fan:
	pn30 = pn3M
	# Reshape into [I1,I2,P]
	P = tail(dim(pn3M),1)
	dim(pn3M) = c(dimb,P)
	dim(data$pns1) = dimb
	dim(data$pns2) = dimb
	# For each fan of each time step, calculate the mean angle of the phase angles, using the magnitude estimates as the lengths of vectors in the direction of the phase angles:
	for(i2 in seq_len(dim(pn3M)[2])){
		for(p in seq_len(dim(pn3M)[3])){
			pn3M[,i2,p] = angleMean(pn3M[,i2,p], data$pns1[,i2] * 10^(data$pns2[,i2]))
			}
		}
	
	# Also, for each fan, calculate the mean angle of the phase angles, using the magnitude estimates as the lengths of vectors in the direction of the phase angles:
	pns3 = NAs(dimb)
	for(i2 in seq_len(dim(pn3M)[2])){
		pns3[,i2] = angleMean(pn3M[,i2,], data$pns1[,i2] * 10^(data$pns2[,i2]))
		}
	dim(pns3) = NULL
	dim(pn3M) = c(prod(dimb),P)
					
	
	
	# Write the mean of the rows of 'out' if required:
	list(pns3=pns3, pn3M=pn3M, pn30=pn30)
}
