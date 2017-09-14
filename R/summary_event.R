#*********************************************
#*********************************************
#' Calculates variaous measures of mean sv, and total backscatter and volume.
#'
#' @param data		A list containing volume and the acoustic data eihter given as 'vbss' for a segment of the data or as 'vbsc'.
#' @param plot.hist	Logical: If TRUE plot the histogram when calculating the robust mean sv.
#' @param minlen	The minimum number of acoustic samples.
#' @param ...		Passed on to density() when calculating the robust mean sv.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom pbapply pblapply
#' @importFrom data.table rbindlist
#'
#' @export
#'
summary_event <- function(event, t=1, plot.hist=FALSE, minlen=1, cores=1, ...){
	summary_one<- function(i, event, plot.hist, minlen){
			vbsc <- read.event(event, var=c("vbsc", "voxels"), t=i) 
			summary_TSD(vbsc, plot.hist=plot.hist, minlen=minlen, ...)
		}
	
	# Get volume equal for all time steps:	
	numt <- read.event(event, var="numt")$numt
	if(identical(t, "all")){
		t <- seq_len(numt)
	}
	
	# Run the summary possibly in parallel:
	availableCores = detectCores()
	if(cores>availableCores){
		warning(paste0("Only ", availableCores, " cores available (", cores, " requested)"))
	}
	cores = min(cores, length(t), availableCores)
	
	if(cores>1){
		cat(paste0("Extracting summary using ", cores, " cores in parallel:\n"))
		# Generate the clusters of time steps:
		cl<-makeCluster(cores)
		# Bootstrap:
		out <- pblapply(t, summary_one, event=event, plot.hist=plot.hist, minlen=minlen, cl=cl)
		# End the parallel bootstrapping:
		stopCluster(cl)
	}
	else{
		cat(paste0("Extracting summary using:\n"))
		out <- pblapply(t, summary_one, event=event, plot.hist=plot.hist, minlen=minlen)
	}
	
	
	data.table::rbindlist(out)
	##################################################
	##################################################
	}



	
		
		
		