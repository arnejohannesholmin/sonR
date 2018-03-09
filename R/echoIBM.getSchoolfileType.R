#*********************************************
#*********************************************
#' Function for extracting the type (dynamic or static) of the school files given by 'schoolfiles'. Four vectors are returned:
#' NA
#'
#' @param schoolfiles  is a vector of paths to school files.
#' @param dynschoolnames  is a vector of four character strings representing dynamic school valiable names.
#' @param staticschoolnames  is a vector of four character strings representing static school valiable names.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSD zeros
#'
#' @export
#' @rdname echoIBM.getSchoolfileType
#'
echoIBM.getSchoolfileType <- function(schoolfiles, dynschoolnames, staticschoolnames, thr=0.5){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-02-20 - Clean version
	########### DESCRIPTION: ###########
	# Function for extracting the type (dynamic or static) of the school files given by 'schoolfiles'. Four vectors are returned:
	#	- schooltypeD, which is 1 for the files that are recognized as dynamic files (containing at least one dynamic variable)
	#	- schooltypeS, which is 1 for the files that are recognized as static files (containing only static variables)
	#	- dynschoolfiles, which are the recognized dynamic files
	#	- staticschoolfiles, which are the recognized static files
	########## DEPENDENCIES: ###########
	# read.TSD()
	############ VARIABLES: ############
	# ---schoolfiles--- is a vector of paths to school files.
	# ---dynschoolnames--- is a vector of four character strings representing dynamic school valiable names.
	# ---staticschoolnames--- is a vector of four character strings representing static school valiable names.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Defining school file type 'schooltype' (0 for dynamic school files and 1 for static school files):
	schooltypeD <- zeros(length(schoolfiles))
	schooltypeS <- zeros(length(schoolfiles))
	schooltypeB <- zeros(length(schoolfiles))
	
	
	########## Execution ##########
	# For loop through the school files:
	for(i in seq_along(schoolfiles)){
		this=suppressWarnings(read.TSD(schoolfiles[i], var=NULL, header=TRUE))
		thislabl <- setdiff(this$labl, c("size","info","d000", labl.TSD("t")))
		# Make an exeption for "size", "info" and "d000", which may be present both in dynamic and static data:
		if(mean(thislabl %in% staticschoolnames) > thr){
			schooltypeS[i] <- 1
			}
		if(mean(thislabl %in% dynschoolnames) > thr){
			schooltypeD[i] <- 1
			}
		if(any(this$labl %in% c("psxS"))){
			schooltypeB[i] <- 1
			}
		}
	
	
	########## Output ##########
	list(schooltypeD=schooltypeD, schooltypeS=schooltypeS, schooltypeB=schooltypeB, dynschoolfiles=schoolfiles[schooltypeD==1], staticschoolfiles=schoolfiles[schooltypeS==1], bothschoolfiles=schoolfiles[schooltypeB==1])
	##################################################
	##################################################
	}
