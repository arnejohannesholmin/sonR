#*********************************************
#*********************************************
#' Estimates high intensity noise in acoustic data by identifying high intensity one dimensional spikes.
#'
#' @param con  is the connection object or a character string naming the output file.
#' @param event  is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
#' @param cruise  is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param t  is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t == "all", all files are read and if t == "none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
#' @param turns  is a numeric giving the length of the runs, causing 'turns' pings to be read at each step in the function.
#' @param k  is a numeric giving the width of the medians across pings, used in runmed() inside get.hins.event_small().
#' @param q  is either a single numeric, or a vector of length 2, givin the start and end point in the linear function along beams defining the threshold times the median filtered values in each direction above which data are classified as high intensity noise.
#' @param beta_school  is a single numeric representing a typical high school value. Only values above 'beta_school' can be classified as high intensity noise. 'beta_school' assures that spikes due to for example fish in very silent regions of the data are not classified as high intensity noise.
#' @param ind  is a list of indexes, as given to subset_TSD(), used to select the subset over which the estimation of high intensity noise is done. Defaulted to exclude the first 100 voxels along each beam.
#' @param max.memory  puts a restiction on the mempry occupied by the function, prompting a call to the user to approve is exceeded.
#' @param ...  variables used in write.TSD() (such as 'ow' for overwriting existing files).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD write.TSD
#' @importFrom utils tail
#'
#' @export
#' @rdname get.hins.event
#'
get.hins.event<-function(con=NULL, event=33, cruise=2009116, esnm="MS70", t=c("20091115084000", "20091115085000"), turns=114, k=5, q=c(1e3, 1e2), beta_school=2e-4, ind=list(-(1:100), NULL), max.memory=2e9, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-06-19 - Clean version.
	# Update: 2010-06-29 - Added the option 'beta_school' and changed to using a two element vector 'q' where the first element is the value of 'q' closest to the sonar, and the second elements is the value of 'q' farthest from the sonar, and other values of 'q' are linearly interpolated between these two values.
	# Update: 2010-07-02 - Changed the number of dimensions of the output 'hini' from 4 to 3.
	# Update: 2010-07-16 - Fixed bug with time steps.
	# Last: 2013-08-21 - Remaned from get.hins.MS70() to get.hins.event().
	########### DESCRIPTION: ###########
	# Estimates high intensity noise in acoustic data by identifying high intensity one dimensional spikes.
	########## DEPENDENCIES: ###########
	# read.event(), zeros(), write.TSD()
	############ VARIABLES: ############
	# ---con--- is the connection object or a character string naming the output file.
	# ---event--- is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
	# ---cruise--- is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
	# ---t--- is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t == "all", all files are read and if t == "none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
	# ---turns--- is a numeric giving the length of the runs, causing 'turns' pings to be read at each step in the function.
	# ---k--- is a numeric giving the width of the medians across pings, used in runmed() inside get.hins.event_small().
	# ---q--- is either a single numeric, or a vector of length 2, givin the start and end point in the linear function along beams defining the threshold times the median filtered values in each direction above which data are classified as high intensity noise.
	# ---beta_school--- is a single numeric representing a typical high school value. Only values above 'beta_school' can be classified as high intensity noise. 'beta_school' assures that spikes due to for example fish in very silent regions of the data are not classified as high intensity noise.
	# ---ind--- is a list of indexes, as given to subset_TSD(), used to select the subset over which the estimation of high intensity noise is done. Defaulted to exclude the first 100 voxels along each beam.
	# ---max.memory--- puts a restiction on the mempry occupied by the function, prompting a call to the user to approve is exceeded.
	# ---...--- variables used in write.TSD() (such as 'ow' for overwriting existing files).
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Function for checking the memory:
	checkMemory<-function(m,M){
		if(m>M){
			ans=readline(paste("Memory exceeding the limit specified (",M,"). Type inn the new limit, or 'y'/'n' to continue/abort",sep=""))
			if(!is.na(as.numeric(ans))){
				M=as.numeric(ans)
				}
			else if(identical(tolower(ans),"y")){
				M=Inf
				}
			else{
				stop("Funciton terminated by the user")
				}
			}
		M
		}
	doubleMemory = 2*.Machine$sizeof.pointer
	
	# Read the time steps of the event:
	time = read.event(event=event, cruise=cruise, esnm="MS70", t=t, var="time")
	if(length(time)==0){
		return(list())
		}
	indt = time$indt
	mtim = time$mtim
	P = length(indt)
	# If the event is empty, write and return empty elements:
	if(P == 0){
		# Define the output variables:
		hins = double()
		hini = array(double(), dim=c(0, 3))
		hinsInd1 = array(double(), dim=c(0, 3))
		hinsInd2 = array(double(), dim=c(0, 3))
		hinsInd3 = array(double(), dim=c(0, 3))
		# Write the data to file:
		if(!is.null(con)){
			write.TSD(list(hins=hins, hini=hini, mtim=mtim), con=con, numt=1, dimension=TRUE, keep.null=TRUE, ...)
			cat("Empty high intensity noise file written for event", event, "\n")
			}
		# Output the data:	
		return(list(hins=hins, hini=hini, indt=indt, mtim=mtim, hinsInd1=hinsInd1, hinsInd2=hinsInd2, hinsInd3=hinsInd3))
		}
	
	if(!all(seq(indt[1], tail(indt, 1)) == indt)){
		stop("The time indexes 't' must correspond to consecutive time steps")
		}
	
	# Read the beam configuration:
	data = read.event(event=event, cruise=cruise, esnm="MS70", var="beams")
	numb = data$numb
	D3 = length(unique(data$freq))
	D2 = numb/D3
	D1 = max(data$lenb)
	
	# Define the runs:
	#runs = ceiling(P/turns)
	#runlengths = c(rep(turns, runs-1), P-(runs-1)*turns)
	#indts = cbind(indt[1]+c(0, cumsum(runlengths[-runs])), indt[1]-1+cumsum(runlengths[]))
	indts = getRuns(turns, indt)
	runs = nrow(indts)
	wideindts = cbind(indts[, 1]-(k-1)/2, indts[, 2]+(k-1)/2)
	wideindts[wideindts<min(indt)] = min(indt)
	wideindts[wideindts>max(indt)] = max(indt)
	
	# Declare the memory count for comparison to 'max.memory':
	memory = 0
	
	# Define the output variables:
	hins = double()
	hini = array(double(), dim=c(0, 3))
	hinsInd1 = array(double(), dim=c(0, 3))
	hinsInd2 = array(double(), dim=c(0, 3))
	hinsInd3 = array(double(), dim=c(0, 3))
	
	# Run through the runs and estimate the high intensity noise:
	cat("Estimating high intensity noise for event ", event, " (", P, " time steps) in ", runs, " runs:\nRun ", sep="")
	for(i in seq_len(runs)){
		cat(i, "\n\n")
		thisindt = indts[i, 1]:indts[i, 2]
		thiswideindt = wideindts[i, 1]:wideindts[i, 2]
		# Locate the high intensity noise for the current run:
		
		thishins = get.hins.event_small(con=con, event=event, cruise=cruise, t=thiswideindt, k=k, q=q, beta_school=beta_school, ind=ind, max.memory=max.memory)
		
		# Check for too large occupied memory:
		memory = memory + sum(sapply(thishins, length))*doubleMemory
		max.memory = checkMemory(memory, max.memory)
		
		# Select only the time steps in 'thisindt', not all used in 'thiswideindt', which include some time steps on either side to comply with the median filtering:
		if(length(thishins$hini)>0){
			valid = thishins$hini[, 3] %in% thisindt
			hini = rbind(hini, thishins$hini[valid, ])
			hins = c(hins, thishins$hins[valid])
			}
		if(length(thishins$hinsInd1)>0){
			valid = thishins$hinsInd1[, 3] %in% thisindt
			hinsInd1 = rbind(hinsInd1, thishins$hinsInd1[valid, ])
			}
		if(length(thishins$hinsInd2)>0){
			valid = thishins$hinsInd2[, 3] %in% thisindt
			hinsInd2 = rbind(hinsInd2, thishins$hinsInd2[valid, ])
			}
		if(length(thishins$hinsInd3)>0){
			valid = thishins$hinsInd3[, 3] %in% thisindt
			hinsInd3 = rbind(hinsInd3, thishins$hinsInd3[valid, ])
			}
		}
	
	
	# Write the data to file:
	if(!is.null(con)){
		write.TSD(list(hins=hins, hini=hini, mtim=mtim), con=con, numt=1, dimension=TRUE, ...)
		cat("High intensity noise file written for event", event)
		}
	
	# Output the data:
	list(hins=hins, hini=hini, indt=indt, mtim=mtim, hinsInd1=hinsInd1, hinsInd2=hinsInd2, hinsInd3=hinsInd3, esnm=esnm)
	##################################################
	##################################################
	}
