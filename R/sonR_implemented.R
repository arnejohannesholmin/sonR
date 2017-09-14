#*********************************************
#*********************************************
#' Test for implemented acosutical systems in echoIBM().
#'
#' @param esnm  is the name of the acoustical instrument, given as a four character string, or a list containing an element named "esnm". May be given in lower case.
#' @param type  is a vector of the types of acoustical instrument to check for implementation in echoIBM. As default the function tests against all implemeneted systems.
#' @param notest  is TRUE to return the implemented systems, in which case no test is made.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD strff
#'
#' @export
#' @rdname sonR_implemented
#'
sonR_implemented<-function(esnm, type=c("MBS","MBE","SBE","OFS"), notest=FALSE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2015-09-09 - Clean version.
	########### DESCRIPTION: ###########
	# Test for implemented acosutical systems in echoIBM().
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---esnm--- is the name of the acoustical instrument, given as a four character string, or a list containing an element named "esnm". May be given in lower case.
	# ---type--- is a vector of the types of acoustical instrument to check for implementation in echoIBM. As default the function tests against all implemeneted systems.
	# ---notest--- is TRUE to return the implemented systems, in which case no test is made.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Here the currently implemented acosutical systems are listed (2015-09-09):
	implemented = list(
		MBS = c("ms70"), 
		MBE = c("me70"), 
		SBE = c("ek60", "ek500"), 
		OFS = c("sx80","sx90","sh80","sh90","su80","su90","sp80","sp90","sn90"))
	types = unlist(implemented[toupper(type)])
	if(notest){
		return(types)
		}
	
	
	##### Execution and output #####
	if(is.list(esnm)){
		esnm = esnm$esnm
		}
	any(strff(types, esnm[1]))
	##################################################
	##################################################
	}