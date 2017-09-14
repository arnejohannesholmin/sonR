#*********************************************
#*********************************************
#' Estimates the mean over time steps of close range noise for each voxel of an underwater acoustic system, using mean() or meanexp.quantile() depending on whether the periodic noise is in effect, and writes the result to file in the TSD format if 'con' is given. The near range noise estimates are smoothed with a kernel smoother given by 'filter'. The variance of the Nadaraya-Watson estimator can be derived straight forward by
#'   var( sum( w_i * x_i ) )
#' = sum( var( w_i * x_i ) )
#' = sum( w_i^2 * var( x_i ) )
#' = sum( w_i^2 * s^2 )
#' = sum( w_i^2 ) * s^2 
#' = sum( (k_i / sum(k_i))^2 ) * s^2 
#' = sum( k_i^2 / sum(k_i)^2 ) * s^2 
#' = sum( k_i^2 ) / sum(k_i)^2 * s^2 
#'
#' @param con  is the connection object or a character string naming the output file.
#' @param event  is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
#' @param cruise  is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param t  is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
#' @param bgns  is either a list containing the vector of background noise 'bgns' of each beam, og the vector itsself.
#' @param filter  is the kernel used when smoothing the near range noise of each beam.
#' @param nsdn  is the number of standard deviation above the mean nrea range noise for each beam, defining the values that are set to 0. The beam is traced from the sonar and out, and at the first value below the threshold, the near range noise is set to zero for all values at and beyond this value.
#' @param farv  is the number of far voxels used for the calclulation of the mean and standard deviation used when assigning zeros to the smoothed near range noise. For farv=13, the 13 farthest non-missing values are used for each beam.
#' @param surface  is a vector of indexes for the beams affected by the surface noise (extrapolated in the near range noise estimation).
#' @param ...  variables used in write.TSD() (such as 'ow' for overwriting existing files).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD NAs write.TSD zeros
#' @importFrom utils tail
#' @importFrom stats dnorm sd
#'
#' @export
#' @rdname get.nrns.event
#'
get.nrns.event<-function(con=NULL,event=NULL,cruise=NULL,t=NULL,bgns=NULL,filter=dnorm(-20:20,sd=4),nsdn=1000,farv=1000,surface=451:500,esnm="MS70",...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-03-12 - Clean version.
	# Update: 2012-05-23 - Changed to estimate the periodic noise as well, using optim() and estimating the background noise and periodic noise simultaneously.
	# Update: 2012-05-23 - Changed to estimate the background noise first, then optimize to estimate the periodic noise.
	# Update: 2012-05-23 - Changed method for estimation of the background noise to using meanexp.quantile() instead of meanexp[dot]conditional().
	# Update: 2012-05-29 - Expanded to several events.
	# Update: 2012-06-13 - Changed method to use the function meanexp.quantile() for beams affected by the periodic noise.
	# Update: 2012-06-13 - Reverted to using simple mean for all beams, since the periodic noise is negligible comared to the extremely high near range noise.
	# Update: 2012-10-26 - Removed defaults for 't', 'event' and 'cruise'.
	# Update: 2012-11-04 - Added the option 'esnm'.
	# Last: 2013-07-19 - Removed 'near.range' and 'ind'.
	########### DESCRIPTION: ###########
	# Estimates the mean over time steps of close range noise for each voxel of an underwater acoustic system, using mean() or meanexp.quantile() depending on whether the periodic noise is in effect, and writes the result to file in the TSD format if 'con' is given. The near range noise estimates are smoothed with a kernel smoother given by 'filter'. The variance of the Nadaraya-Watson estimator can be derived straight forward by
	#   var( sum( w_i * x_i ) )
	# = sum( var( w_i * x_i ) )
	# = sum( w_i^2 * var( x_i ) )
	# = sum( w_i^2 * s^2 )
	# = sum( w_i^2 ) * s^2 
	# = sum( (k_i / sum(k_i))^2 ) * s^2 
	# = sum( k_i^2 / sum(k_i)^2 ) * s^2 
	# = sum( k_i^2 ) / sum(k_i)^2 * s^2 
	########## DEPENDENCIES: ###########
	# read.event(), zeros(), write.TSD()
	############ VARIABLES: ############
	# ---con--- is the connection object or a character string naming the output file.
	# ---event--- is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
	# ---cruise--- is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
	# ---t--- is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
	# ---bgns--- is either a list containing the vector of background noise 'bgns' of each beam, og the vector itsself.
	# ---filter--- is the kernel used when smoothing the near range noise of each beam.
	# ---nsdn--- is the number of standard deviation above the mean nrea range noise for each beam, defining the values that are set to 0. The beam is traced from the sonar and out, and at the first value below the threshold, the near range noise is set to zero for all values at and beyond this value.
	# ---farv--- is the number of far voxels used for the calclulation of the mean and standard deviation used when assigning zeros to the smoothed near range noise. For farv=13, the 13 farthest non-missing values are used for each beam.
	# ---surface--- is a vector of indexes for the beams affected by the surface noise (extrapolated in the near range noise estimation).
	# ---...--- variables used in write.TSD() (such as 'ow' for overwriting existing files).
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(is.list(bgns)){
		bgns = bgns$bgns
		}
	if(length(bgns)==0){
		warning("Background noise, independent on range (long range), is not given and will be set to 0")
		bgns = 0
		}
	
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
		cruise = rep(cruise,length.out=nevents)
		}
	
	# Read the length of the sample intervals, the number of beams, and the length of the beams:
	sint = NULL
	numb = NULL
	lenb = NULL
	# Read the time steps of the events:
	indt = vector("list",nevents)
	utim = NULL
	for(ev in seq_len(nevents)){
		data = read.event(event=event[ev],cruise=cruise[ev],var="beams",esnm=esnm)
		sint = c(sint,data$sint)
		numb = c(numb,data$numb)
		lenb = c(lenb,max(data$lenb))
		time = read.event(event=event[ev],cruise=cruise[ev],t=t[[ev]],var="time",esnm=esnm)
		indt[[ev]] = time$indt
		utim = c(utim,time$utim)
		}
	if(length(unlist(indt))==0){
		stop("No time steps chosen")
		}
	ntimesteps = length(unlist(indt))
	
	# Check if the beam configuration is equal across events:
	if(!all(sint==sint[1]) || !all(numb==numb[1]) || !all(lenb==lenb[1])){
		stop("Events must have the same beam configuration")
		}
	else{
		sint = sint[1]
		numb = numb[1]
		lenb = lenb[1]
		}
	
	# If 'bgns' has length shorter than the number of beams, repeat to the required length:
	if(length(bgns)<numb){
		bgns = rep(bgns,length.out=numb)
		}
	
	# Define the farthest voxel, used for extracting the mean and standard deviation of the near range noise. From these a limit value of the smoothed near range noise is defined, by the mean + nsdn * sd * sqrt(sum(filter^2)/sum(filter)^2). This can be derived by taking the variance of the Nadaraya-Watson estimator (see desciption at the top).
	if(farv<1){
		farv = floor(lenb*farv)
		}
	else if(farv>=lenb){
		stop("The specifier for the far voxels used in the calculation of the mean and standard deviation of the near range noise, in turn used for defining values as 0, beyond the first value below the 0-threshold (mean + nsdn * sd * sum(filter^2)/sum(filter)^2), is larger than the number of voxels along beams")
		}
		
	
	##### Execution and output #####
	### 1. Estimate the near range noise for each voxel by the simple mean over all time steps for the healthy beams and by the mean of the lower quantile for beams affected by the periodic noise: ###
	nrns = zeros(lenb,numb)
	
	# Run through the events, and through the time steps in the events:
	for(ev in seq_len(nevents)){
		cat(paste("Estimating near range noise, event ",event[ev],"\nPings:\n"))
		for(i in indt[[ev]]){
			cat(i," ",sep="")
			# Import acoustic data:
			Sv = read.event(event=event[ev],cruise=cruise[ev],t=i,var="vbsc",TVG=FALSE,esnm=esnm)$vbsc
			# Use the simple mean:
			nrns[seq_len(nrow(Sv)),] = nrns[seq_len(nrow(Sv)),] + Sv
			}
		cat("\n")
		}
	# Divide by the number of time steps to obtain the mean, and subtract the background noise:
	nrns = nrns/ntimesteps - matrix(bgns,nrow=lenb,ncol=numb,byrow=TRUE)
	
	# Smooth the near range noise and assign 0 to values outside a range specified by the mean and standard deviation of the filtered version under the hypothesis of zero trend. (See description at the top):
	sdfact = sqrt(sum(filter^2)/sum(filter)^2)
	# Declare the smoothed near range noise, and the smoothed near range noise inserted zeros at and beyond the first occurance of a value below the threshold.
	nrn0 = nrns
	snrn = nrns
	snr0 = nrns
	thrn = NAs(numb)
	thrs = NAs(numb)
	zern = NAs(numb)
	zers = NAs(numb)
	
	for(i in seq_len(numb)){
		# Select only the non-missing values:
		notNA = !is.na(nrns[,i])
		# Apply the smoothing filter to the near range noise estimates:
		snrn[notNA,i] = filter2(nrns[notNA,i],filter)
		snr0[notNA,i] = snrn[notNA,i]
		# Get the mean and standard deviation of the values beyond some voxel specified by 'farv':
		meannrns = mean(tail(nrns[notNA,i],farv),na.rm=TRUE)
		sdnrns = sd(tail(nrns[notNA,i],farv),na.rm=TRUE)
		sdsnrn = sdnrns*sdfact
		thrs[i] = meannrns+sdsnrn*nsdn
		# Insert zeros:
		zers[i] = max(which(snrn[,i]>thrs[i]))
		snr0[seq(zers[i],lenb),i] = 0
		
		# Add zeros in a similar way to the unsmoothed near range noise:
		thrn[i] = meannrns+sdnrns*nsdn
		zern[i] = max(which(nrns[,i]>thrn[i]))
		nrn0[seq(zern[i],lenb),i] = 0
		}
	
	# If 'surface' is given, extrapolate to estimate the surface beams:
	if(length(surface)>0 && all(surface)>=0){
		if(all(c("dira","dire") %in% names(data))){
			for(i in surface){
				ind_i = which.min(sqrt((data$dira[-surface]-data$dira[i])^2+(data$dire[-surface]-data$dire[i])^2))
				nrns[,i] = nrns[,ind_i]
				nrn0[,i] = nrn0[,ind_i]
				snrn[,i] = snrn[,ind_i]
				snr0[,i] = snr0[,ind_i]
				}
			}
		#else if(all(c("dirl","dirt") %in% names(data))){
		#	for(i in surface){
		#		ind_i = which.min(sqrt((data$dirl[-surface]-data$dirl[i])^2+(data$dirt[-surface]-data$dirt[i])^2))
		#		nrns[,i] = nrns[,ind_i]
		#		nrn0[,i] = nrn0[,ind_i]
		#		snrn[,i] = snrn[,ind_i]
		#		snr0[,i] = snr0[,ind_i]
		#		}
		#	}
		else{
			#warning("Could not extrapolate to the specified surface beams. Neither \"dira\",\"dire\" nor \"dirl\",\"dirt\" given")
			warning("Could not extrapolate to the specified surface beams (\"dira\" or \"dire\" not given)")
			}
		}
	
	# Write the mean of the rows of 'out' if required:
	out = list(nrns=nrns,nrn0=nrn0,snrn=snrn,snr0=snr0,fltr=filter,thrn=thrn,thrs=thrs,zern=zern,zers=zers,utim=utim,nsdn=nsdn,sfni=surface,esnm=esnm)
	if(!is.null(con)){
		write.TSD(out,con=con,numt=ntimesteps,...)
		}
	out
	##################################################
	##################################################
	}
