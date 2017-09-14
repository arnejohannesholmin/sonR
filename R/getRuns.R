#*********************************************
#*********************************************
#' Generates runs of length 'n' of the vector 't' of time steps.
#'
#' @param n  is the length of the runs.
#' @param t  is a vector of time steps to group into runs.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname getRuns
#'
getRuns<-function(n,t){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-08-21 - Clean version.
	########### DESCRIPTION: ###########
	# Generates runs of length 'n' of the vector 't' of time steps.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---n--- is the length of the runs.
	# ---t--- is a vector of time steps to group into runs.
	
	
	##################################################
	##################################################
	# Get the number of runs:
	T=length(t)
	runs=ceiling(T/n)
	# Get the run lengths, constant up to the last run:
	runlengths=c(rep(n,runs-1),T-(runs-1)*n)
	# Generate the runs from the 'indt' input:
	cbind(t[1]+c(0,cumsum(runlengths[-runs])), t[1]-1+cumsum(runlengths[]))
	##################################################
	##################################################
	}

