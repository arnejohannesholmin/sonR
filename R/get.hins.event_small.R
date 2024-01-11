#*********************************************
#*********************************************
#' Estimates high intensity noise in an underwater acoustic system, by identifying high intensity one dimensional spikes.
#'
#' @param con  is the connection object or a character string naming the output file.
#' @param event  is the identifier of the event. Given either as the number of the event, a string contained in the name of the event, or the path of the event directory. More than one event may be given, in a vector.
#' @param cruise  is the identifier of the cruise. Given either as the specification used by IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param t  is the identifier of the time points. Given either as a vector of integers between 1 and the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points are extracted. If t=="all", all files are read and if t=="none" no action is done. If more than one event is given, 't_bgns' must be given as a list of vectors.
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
#' @importFrom TSD arr.ind2ind ind.expand ind2arr.ind NAs zeros
#' @importFrom stats runmed
#'
#' @export
#' @rdname get.hins.event_small
#'
get.hins.event_small<-function(con=NULL, event=1, cruise=2009116, esnm="MS70", t=1:10, k=5, q=c(1e3, 1e2), beta_school=2e-4, ind=list(-(1:100), NULL), max.memory=2e9, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-06-19 - Clean version.
	# Update: 2010-06-29 - Added the option 'beta_school' and changed to using a two element vector 'q' where the first element is the value of 'q' closest to the sonar, and the second elements is the value of 'q' farthest from the sonar, and other values of 'q' are linearly interpolated between these two values.
	# Update: 2010-07-02 - Changed the number of dimensions of the output 'hini' from 4 to 3.
	# Update: 2010-07-10 - Fixed bug with displaced time steps reported, and problem with Sv not being segmented with 'ind' during extractino of means.
	# Update: 2010-07-16 - Fixed bug with time steps.
	# Last: 2013-08-21 - Remaned from get.hins.MS70_small() to get.hins.event_small().
	
	
	##################################################
	##################################################
	##### Preparation #####
	doubleMemory = 2*.Machine$sizeof.pointer
	
	# Read the time steps of the event:
	time = read.event(event=event, cruise=cruise, esnm=esnm, t=t, var="time")
	indt = time$indt
	mtim = time$mtim
	P = length(indt)
	# Read the beam configuration:
	data = read.event(event=event, cruise=cruise, esnm=esnm, var="beams", t=indt[1])
	numb = data$numb[1]
	D3 = length(unique(c(data$freq)))
	D2 = numb/D3
	D1 = max(data$lenb)
	# Define the q-sequence, which is a linear function of the range, and the mean q, used for the along beam means:
	q = rep(q, length.out=2)
	qseq = seq(q[1], q[2], length.out=D1)
	
	# Get the indexes specified by 'ind' as typed into []:
	ind1 = ind.expand(ind[c(2, 3)], c(D2, D3))
	ind2 = ind.expand(ind[c(1, 3)], c(D1, D3))
	ind3 = ind.expand(ind[c(1, 2)], c(D1, D2))
	
		
	##### Execution and output #####
	# Error if too much memory is required to run the function:
	total.memory = (prod(D2, D3, P) + prod(D1, D3, P) + prod(D1, D2, P)) * 2*doubleMemory
	if((total.memory) > max.memory){
		stop(paste("Too much memory required (", total.memory/1e6, " MB) . Try reducing the number of time steps or increasing the value of 'max.memeory'", sep=""))
		}
	
	# Define the arrays holding the means along each dimension:
	dimMeans1 = zeros(D2, D3, P)
	dimMeans2 = zeros(D1, D3, P)
	dimMeans3 = zeros(D1, D2, P)
	
	# Extract means along each dimension:
	cat("Extracting means along each dimension:\nPings:\n")
	for(i in seq_len(P)){
		cat(indt[i], " ", sep="")
		# Import acoustic data and ignore near-range and apply the subset specified by 'ind':
		Sv = read.event(event=event, cruise=cruise, esnm=esnm, t=indt[i], var=c("vbsc", "lenb"), TVG=TRUE)
		tempD1 = max(Sv$lenb)
		tempind2 = ind.expand(ind[c(1, 3)], c(tempD1, D3))
		
		Sv = Sv$vbsc
		dim(Sv) = c(tempD1, D2, D3)
		# Get the means along the dimensions, appliyng the subset along beams (but not across beams or time steps):
		#dimMeans1[, , i]=apply(Sv[tempind2[[1]], , , drop=FALSE], c(2, 3), mean, na.rm=TRUE)
		#dimMeans2[, , i]=apply(Sv, c(1, 3), mean, na.rm=TRUE)
		#dimMeans3[, , i]=apply(Sv, c(1, 2), mean, na.rm=TRUE)
		
		temp = Sv[tempind2[[1]], , , drop=FALSE]
		dim(temp) = c(dim(temp)[1], prod(dim(temp)[2:3]))
		dimMeans1[, , i] = colMeans(temp)
		
		temp = aperm(Sv,c(2,1,3))
		dim(temp) = c(dim(temp)[1], prod(dim(temp)[-1]))
		dimMeans2[seq_len(tempD1), , i] = colMeans(temp)
		
		temp = Sv
		dim(temp) = c(prod(dim(temp)[1:2]), dim(temp)[3])
		dimMeans3[seq_len(tempD1), , i] = rowMeans(temp)
		}
	cat("\n")
		
	# Identify and discard the missing values:
	NA1 = apply(dimMeans1, 1:2, function(x) any(is.na(x)))
	NA2 = apply(dimMeans2, 1:2, function(x) any(is.na(x)))
	NA3 = apply(dimMeans3, 1:2, function(x) any(is.na(x)))
	
	# Get median across pings:
	medians1 = NAs(D2, D3, P)
	medians2 = NAs(D1, D3, P)
	medians3 = NAs(D1, D2, P)
	cat("Running median filter along first dimension...\n")
	for(i2 in seq_len(D2)){
		for(i3 in seq_len(D3)){
			if(!NA1[i2, i3]){
				medians1[i2, i3, ] = runmed(dimMeans1[i2, i3, ], k)
				}
			}
		}
	cat("Running median filter along second dimension...\n")
	for(i1 in seq_len(D1)){
		for(i3 in seq_len(D3)){
			if(!NA2[i1, i3]){
				medians2[i1, i3, ] = runmed(dimMeans2[i1, i3, ], k)
				}
			}
		}
	cat("Running median filter along third dimension...\n")
	for(i1 in seq_len(D1)){
		for(i2 in seq_len(D2)){
			if(!NA3[i1, i2]){
				medians3[i1, i2, ] = runmed(dimMeans3[i1, i2, ], k)
				}
			}
		}
	
	# Identify the high intensity noise as the voxels that exceed 'q' times the median of the k neighboring time steps for that voxel:
	hinsInd1 = which(dimMeans1/medians1>max(q) & dimMeans1>beta_school , arr.ind=TRUE)
	hinsInd2 = which(dimMeans2/medians2>qseq & dimMeans2>beta_school , arr.ind=TRUE)
	hinsInd3 = which(dimMeans3/medians3>qseq & dimMeans3>beta_school , arr.ind=TRUE)
	# Discard the voxels not included in 'ind':
	hinsInd1 = hinsInd1[hinsInd1[, 1] %in% ind1[[1]] & hinsInd1[, 2] %in% ind1[[2]], , drop=FALSE]
	hinsInd2 = hinsInd2[hinsInd2[, 1] %in% ind2[[1]] & hinsInd2[, 2] %in% ind2[[2]], , drop=FALSE]
	hinsInd3 = hinsInd3[hinsInd3[, 1] %in% ind3[[1]] & hinsInd3[, 2] %in% ind3[[2]], , drop=FALSE]
	# Add the first time steps to the third column to obtain actual time steps (not only time step indexes within 't').
	hinsInd1[, 3] = hinsInd1[, 3]+t[1]-1
	hinsInd2[, 3] = hinsInd2[, 3]+t[1]-1
	hinsInd3[, 3] = hinsInd3[, 3]+t[1]-1
	# Discard the along beams high intensity noise if across beam high intensity noise is detected:
	if(length(hinsInd2)>0 | length(hinsInd3)>0){
		hinsInd1 = zeros(0, 3)
		}
	
	# Warning if the value of 'q' is too low, which will require too much memory:
	total.memory = (nrow(hinsInd1)*D1 + nrow(hinsInd2)*D2 + nrow(hinsInd3)*D3) *3*doubleMemory
	if((total.memory) > max.memory){
		warning("Voxel indexes not returned")
		list(indt=indt, mtim=mtim, hinsInd1=hinsInd1, hinsInd2=hinsInd2, hinsInd3=hinsInd3)
		}
	# Expand to voxel indexes:
	else{
		# Uniqueify the voxels that are identified as high intensity noise (voxels may be included in beams along several dimensions)
		c11 = rep(seq_len(D1), each=nrow(hinsInd1))
		c21 = hinsInd1[, 1]
		c31 = hinsInd1[, 2]
		c41 = hinsInd1[, 3]
		c12 = hinsInd2[, 1]
		c22 = rep(seq_len(D2), each=nrow(hinsInd2))
		c32 = hinsInd2[, 2]
		c42 = hinsInd2[, 3]
		c13 = hinsInd3[, 1]
		c23 = hinsInd3[, 2]
		c33 = rep(seq_len(D3), each=nrow(hinsInd3))
		c43 = hinsInd3[, 3]
		uniquehinsInd = unique(  c( arr.ind2ind(cbind(c11, c21, c31, c41), c(D1, D2, D3, max(indt))),   arr.ind2ind(cbind(c12, c22, c32, c42), c(D1, D2, D3, max(indt))),   arr.ind2ind(cbind(c13, c23, c33, c43), c(D1, D2, D3, max(indt))) )  )
		# Convert to array indexes:
		uniquehinsInd = ind2arr.ind(uniquehinsInd, c(D1, D2*D3, max(indt)))
		hinsP = sort(unique(uniquehinsInd[, 3]))
		
		# Define the vector holding the high intensity noise values:
		uniquehins = NAs(nrow(uniquehinsInd))
		
		# Extract the high intensity noise at the identified voxels by moving through the time steps:
		cat("Extracting high intensity noise values:\nPings:\n")
		for(i in seq_along(hinsP)){
			cat(hinsP[i], " ", sep="")
			# Import acoustic data and ignore near-range and apply the subset specified by 'ind':
			Sv = read.event(event=event, cruise=cruise, esnm=esnm, t=hinsP[i], var="vbsc", TVG=FALSE)$vbsc
			# Get the high intensity noise:
			this = which(uniquehinsInd[, 3] %in% hinsP[i])
			uniquehins[this] = Sv[arr.ind2ind(uniquehinsInd[this, 1:2], c(D1, D2*D3))]
			}
		cat("\n")
		# Remove the missing values (high intensity noise should not be estimated at the missing values):
		nas = is.na(uniquehins)
		if(any(nas)){
			uniquehins = uniquehins[!nas]
			uniquehinsInd = uniquehinsInd[!nas, ]
			}
		
		list(hins=uniquehins, hini=uniquehinsInd, indt=indt, mtim=mtim, hinsInd1=hinsInd1, hinsInd2=hinsInd2, hinsInd3=hinsInd3, esnm=esnm)
		}
	##################################################
	##################################################
	}
