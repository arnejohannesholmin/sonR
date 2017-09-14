#*********************************************
#*********************************************
#' Removes duplicated time steps from the list of time step indices 'tlist'.
#'
#' @param tlist  is a list of time steps for each file in the vector of file paths 'files'.
#' @param files  is a list of file names.
#' @param ind  is a vector of indices for the files for which to check for duplicated time steps.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname removeDuplicated_tlist
#'
removeDuplicated_tlist<-function(tlist, files, ind=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-10-01 - Clean version.
	########### DESCRIPTION: ###########
	# Removes duplicated time steps from the list of time step indices 'tlist'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---tlist--- is a list of time steps for each file in the vector of file paths 'files'.
	# ---files--- is a list of file names.
	# ---ind--- is a vector of indices for the files for which to check for duplicated time steps.
	

	##################################################
	##################################################
	##### Preparation #####
	if(length(ind)==0){
		ind = seq_along(tlist)
		}
	# The list 'duptlist' contains logical values for each time step of each file to be read, identifying duplicated time steps for the school variables:
	duptlist = split(duplicated(unlist(tlist[ind], use.names=FALSE)), rep(seq_along(tlist[ind]), sapply(tlist[ind], length)))
	
	
	##### Execution #####
	if(any(unlist(duptlist, use.names=FALSE)) && all(sapply(duptlist,length)>1)){
		warning(paste("Duplicated time steps found and removed from the output:", paste(files[ind][sapply(duptlist,sum)>0], collapse="\n"), sep="\n"))
		for(i in seq_along(duptlist)){
			if(length(tlist[ind][[i]])>0){
				tlist[ind][[i]] = setdiff(tlist[ind][[i]], which(duptlist[[i]]))
				}
			}
		}
	
	
	##### Output #####
	tlist
	##################################################
	##################################################
	}
