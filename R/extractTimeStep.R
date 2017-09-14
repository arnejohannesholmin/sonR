#*********************************************
#*********************************************
#' Extracts a time step from TSD-data, so that if a one dimensional variable has two dimensions, the second is regarded as time and the apropriate time steps i extracted.
#'
#' @param data  is the list of TSD inputs as returned from read.TSD.
#' @param t  is a single integer giving the time step to extract.
#' @param var  is a vector of the TSD-names og the variable to extract the time step from. Variables that are not present in 'var' are returned unaltered. The default var='all' select all variables.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD labl.TSD
#'
#' @export
#' @rdname extractTimeStep
#'
extractTimeStep<-function(data, t=1, var="all"){
		
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2015-04-27 - Clean version.
	########### DESCRIPTION: ###########
	# Extracts a time step from TSD-data, so that if a one dimensional variable has two dimensions, the second is regarded as time and the apropriate time steps i extracted.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is the list of TSD inputs as returned from read.TSD.
	# ---t--- is a single integer giving the time step to extract.
	# ---var--- is a vector of the TSD-names og the variable to extract the time step from. Variables that are not present in 'var' are returned unaltered. The default var='all' select all variables.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Get indices of the elements to treat
	if(identical(var,"all")){
		var = seq_along(data)
		}
	else{
		var = which(names(data) %in% var)
		}
	
	
	##### Execution #####
	# Only consider the veriables with positive length:
	emptyvar = intersect(var, which(sapply(data, length)==0))
	var = setdiff(var, emptyvar)
	# Get the list variables, and extract the time step from the lists here:
	listvar = intersect(var, which(sapply(data, is.list)))
	data[listvar] = lapply(data[listvar], function(y) y[[t]])
	# Update 'var' to exclude the list variables:
	var = setdiff(var, listvar)
	
	# Define the variables expected to have one, two or three dimensions when time is along the last dimension:
	var1 = intersect(var, which(names(data) %in% c(labl.TSD(var="v"), labl.TSD(var="t"), labl.TSD(var="b1"))))
	var2 = intersect(var, which(names(data) %in% labl.TSD(var="b2")))
	var3 = intersect(var, which(names(data) %in% c("vbsc","mvbs","angl","angt","volx","harx","psxx","psyx","pszx")))
	
	# Get the number of time steps for each variable
	var1 = var1[sapply(data[var1], length) > 1]
	var2 = var2[sapply(data[var2], function(y) length(dim(y))) == 2]
	var3 = var3[sapply(data[var3], function(y) length(dim(y))) == 3]
	
	# Extract time steps:
	data[var1] = lapply(data[var1], function(y) y[t])
	data[var2] = lapply(data[var2], function(y) y[,t])
	data[var3] = lapply(data[var3], function(y) y[,,t])
	
	
	##### Output #####
	# Return a list of the edgepoints of the circular voxels:
	data
	##################################################
	##################################################
	}
