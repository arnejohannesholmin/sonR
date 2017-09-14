#*********************************************
#*********************************************
#' Shifts the vessel information 'rtzv', 'lonv', 'latv', and 'ispv' by the time offset 'TOV'. Used in read.event().
#'
#' @param var  is a vector of the names of the list holding the vessel data (obtained by names(out) in read.event()).
#' @param time  is a list of the time information (obtained by out[c("mtim","utim")] in read.event()).
#' @param rawvessel  is a vector of the paths to the TSD files from which the data are read.
#' @param TOV  is the time offset of the vessel information, caused by some error in the SIMRAD system.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ang2rot mtim.TSD read.TSDs utim.TSD
#' @importFrom stats approx
#'
#' @export
#' @rdname correctVessel
#'
correctVessel<-function(var, time, rawvessel, TOV=0){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-05-06 - Clean version.
	########### DESCRIPTION: ###########
	# Shifts the vessel information 'rtzv', 'lonv', 'latv', and 'ispv' by the time offset 'TOV'. Used in read.event().
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---var--- is a vector of the names of the list holding the vessel data (obtained by names(out) in read.event()).
	# ---time--- is a list of the time information (obtained by out[c("mtim","utim")] in read.event()).
	# ---rawvessel--- is a vector of the paths to the TSD files from which the data are read.
	# ---TOV--- is the time offset of the vessel information, caused by some error in the SIMRAD system.
	
	
	##################################################
	##################################################
	##### Preparation #####
	rawvariables=c("rtzv","lonv","latv","ispv")
	irawvariables=c("irzv","ilnv","iltv","iisv")
	rawinvar=intersect(rawvariables,var)
	thisout=list()
	
	
	##### Execution and output #####
	# If any of 'rawvariables' are requested, read the time shifted information from the raw dynamic vessel information file:
	if(length(rawinvar)>0){
		ipresentnames = irawvariables[rawinvar %in% rawvariables]
		presentnames = rawvariables[rawinvar %in% rawvariables]
		# Read the raw vessel file:
		if(!is.list(rawvessel)){
			rawvessel = read.TSDs(rawvessel, var=c("iutm","imtm",ipresentnames), t="all")
			}
		
		# If rotation angles are given, convert to rotation angles:
		if(any(ipresentnames == "irzv")){
			rawvessel$irzv = ang2rot(rawvessel$irzv)
			}
			
		# Locate the time shifted vessel information:
		if(length(rawvessel$iutm)>0){
			for(i in seq_along(ipresentnames)){
				thisout[[presentnames[i]]] = approx(x=rawvessel$iutm, y=rawvessel[[ipresentnames[i]]], xout=unique(unlist(utim.TSD(time[c("mtim","utim")])))+TOV, method="linear", rule=2)$y
				}
			}
		else if(length(rawvessel$imtm)>0){
			for(i in seq_along(ipresentnames)){
				thisout[[presentnames[i]]] = approx(x=rawvessel$imtm, y=rawvessel[[ipresentnames[i]]], xout=unique(unlist(mtim.TSD(time[c("mtim","utim")])))+TOV/86400, method="linear", rule=2)$y
				}
			}
		
		# Convert back to angles in the range [0, 2*pi]:
		if(any(ipresentnames == "irzv")){
			thisout$rtzv = thisout$rtzv %% (2*pi)
			}
		}
	thisout
	##################################################
	##################################################
	}
