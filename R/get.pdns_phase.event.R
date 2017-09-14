#*********************************************
#*********************************************
#' Estimates the phase of the periodic noise of the MS70 sonar. If the variables 'badb', 'pns1', 'pns2', 'harm' and 'bgns' are not present in the input list 'data', an attempt is made to read them from the directory holding the MS70 noise estimate for the given event. Only used in echoIBM.vbsc2p.event().
#'
#' @param event  is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
#' @param cruise  is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param t  is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
#' @param noise  is a list of the periodic noise parameters for the magnitude (pns1), kurtosis (pns2), harmonics (harm) and identifyers for the periodic beams (badb), and the background noise (bgns) for use in the estimation of periodic noise. If not given, these are attempted to be read from the noise files of the cruise.
#' @param nsind  is a vector of indexes along the beams, as input to ind.expand(), used to select the subset over which the estimation of the phase of the periodic noise is done. If given as a single numeric, the outermost 'nsind' voxels are used in each beam.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param pdns_scale  is used to scale the noise in order to allow the optimization to work.
#' @param TVG  is FALSE if TVG compensation is to be removed from the data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname get.pdns_phase.event
#'
get.pdns_phase.event<-function(event=1, cruise=2009116, esnm="MS70", t=1:100, noise=list(), nsind=0.75, dir.data=NULL, pdns_scale=1e-14, TVG=TRUE){
	
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
	# Update: 2012-06-13 - Changed to identify the beams affected by the periodic noise by storing the relative height of the highest peak in spectral density compared to the specral density of the zero-frequency, averaging over all pings, and thresholding with 'fftthreshold'.
	# Update: 2012-06-30 - Changed the treatment of 'ind1' and 'near.range' to using ind.expand().
	# Update: 2012-06-30 - Changed to use magnitude*10^kurtosis as weights for the estimation of the phase.
	# Update: 2013-08-23 - Changed to use get.pdns_phase.TSD().
	# Last: 2013-07-19 - Changed names from get.pdns_phase.MS70() to get.pdns_phase.event().
	########### DESCRIPTION: ###########
	# Estimates the phase of the periodic noise of the MS70 sonar. If the variables 'badb', 'pns1', 'pns2', 'harm' and 'bgns' are not present in the input list 'data', an attempt is made to read them from the directory holding the MS70 noise estimate for the given event. Only used in echoIBM.vbsc2p.event().
	########## DEPENDENCIES: ###########
	# read.event(), get.pdns_phase.TSD()
	############ VARIABLES: ############
	# ---event--- is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
	# ---cruise--- is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
	# ---t--- is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
	# ---noise--- is a list of the periodic noise parameters for the magnitude (pns1), kurtosis (pns2), harmonics (harm) and identifyers for the periodic beams (badb), and the background noise (bgns) for use in the estimation of periodic noise. If not given, these are attempted to be read from the noise files of the cruise.
	# ---nsind--- is a vector of indexes along the beams, as input to ind.expand(), used to select the subset over which the estimation of the phase of the periodic noise is done. If given as a single numeric, the outermost 'nsind' voxels are used in each beam.
	# ---dir.data--- is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
	# ---pdns_scale--- is used to scale the noise in order to allow the optimization to work.
	# ---TVG--- is FALSE if TVG compensation is to be removed from the data.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Read the beams variables:
	data = c(noise,read.event(event=event,cruise=cruise,esnm=esnm,var="beams"))
	# Add the time of the event:
	time=read.event(event=event,cruise=cruise,esnm=esnm,t=t,var="time")
	data[names(time)]=time
	
	
	##### Execution and output #####
	get.pdns_phase.TSD(data, nsind=nsind, cruise=cruise, event=event, esnm=esnm, dir.data=NULL, pdns_scale=pdns_scale, TVG=TVG)
	##################################################
	##################################################
	}
