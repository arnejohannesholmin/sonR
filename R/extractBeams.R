#*********************************************
#*********************************************
#' Extracts a beams from TSD-data, so that if a one dimensional variable has two dimensions, the first is regarded as beams and the apropriate time beams are extracted.
#'
#' @param data  is the list of TSD inputs as returned from read.TSD.
#' @param beams  is an integer vectir giving the beams to extract.
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
#' @rdname extractBeams
#'
extractBeams <- function(data, beams=1, var="all", drop=FALSE){
		
	############### LOG: ###############
	# Start: 2018-03-09 - Modified extractTimeStep()
	
	
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
	data[listvar] = lapply(data[listvar], function(y) y[[beams]])
	# Update 'var' to exclude the list variables:
	var = setdiff(var, listvar)
	
	# Define the variables expected to have one, two or three dimensions when time is along the last dimension:
	var2 = intersect(var, which(names(data) %in% labl.TSD(var="b2")))
	var3 = intersect(var, which(names(data) %in% c("vbsc","mvbs","angl","angt","volx","harx","psxx","psyx","pszx")))
	
	# Get the number of time steps for each variable
	var2 = var2[sapply(data[var2], function(y) length(dim(y))) == 2 & sapply(data[var2], NROW) >= max(beams)]
	var3 = var3[sapply(data[var3], function(y) length(dim(y))) == 3 & sapply(data[var3], function(y) dim(y)[2]) >= max(beams)]
	
	# Extract time steps:
	data[var2] = lapply(data[var2], function(y) y[beams, , drop=drop])
	data[var3] = lapply(data[var3], function(y) y[, max(beams), , drop=drop])
	
	
	##### Output #####
	data
	##################################################
	##################################################
	}
