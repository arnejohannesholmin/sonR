#*********************************************
#*********************************************
#' Locates the paths to files holding the noise estimates and the estimates of the correlation between voxels. First the relevant files are seached for in the event given by 'event', 'cruise', and 'esnm', and then the files are seached for in the echoIBM/Resources/Noise/Main directory by matching 'esnm' and 'utim', or by 'esnm', 'cruise' if 'utim' is not given. This may select more than one relevant file, and read.TSDs() should be applied to the detected files, possibly with clean=FALSE.
#'
#' @param cruise  is either the idenfication number of the cruise (or to the directory "/Applications/echoIBM/Resources" in which default noise files for echoIBM are located), given as specified by the IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param utim  is the unix time by which the noise files should be selected.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param noisevar  is a vetor of TSD names which mest be present in the noise estimate files.
#' @param corvar  is a vetor of TSD names which mest be present in the correlation estimate files.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname noise.path
#'
noise.path<-function(event=NULL, cruise=2009116, esnm="MS70", utim=NULL, dir.data=NULL, noisevar=c("bgns","nrnp"), corvar=c("crb1","olpb")){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-07-07 - Clean version.
	# Update: 2012-09-25 - Removed the selection of the first file, if more than one file is present. As a consequence, read.TSDs() should be applied on the output instead of read.TSD().
	# Last: 2012-12-06 -Changed to first look for the noise in the event, and then in the directory "/Applications/echoIBM/Resources".
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# If not given, set the data directory containing data (structured in the directory structure specified in the documentation of echoIBM) as the string Acoustics_datasets_directory().
	if(is.null(dir.data)){
		dir.data=Acoustics_datasets_directory()
		}
	
	# Get the path to the event
	if(length(event)>0){
		event = event.path(event=event,cruise=cruise,dir.data=dir.data)$event
		cruisedir = file.path(dirname(dirname(dirname(dirname(event)))), "Noise", "Main")
		}
	else{
		cruisedir = NULL
	}
	# Define the location of the commonly stored files:
	#commondir=file.path(Acoustics_datasets_directory(),"Resources","Noise","Main")
	
	
	########## Execution ##########
	out=list()
	# (1) Try searching for the desired variables in the event:
	if(length(event)>0){
		out=locate.noise_file(event,noisevar,corvar,esnm)
	}
	
	if(length(out)==0 && length(cruisedir)>0){
		out=locate.noise_file(cruisedir,noisevar,corvar,esnm)
	}

	if(length(out)==0){
		stop("Global noise files no longer supported. Noise files should be put in the Event directory.")
		#out=locate.noise_file(commondir,noisevar,corvar,esnm)
	}

	# Try searching for the noise in the directory file.path(dirname(echoIBM_frameworks),"Resources","Noise","Main"):
	if(any(length(out$noisefiles)+length(out$corfiles)==0)){
		#thisout=locate.noise_file(file.path(dirname(echoIBM_frameworks),"Resources","Noise","Main"),noisevar,corvar,esnm)
		thisout=locate.noise_file(cruisedir,noisevar,corvar,esnm)
		if(length(thisout)==0){
			#thisout=locate.noise_file(commondir,noisevar,corvar,esnm)
			}
		
		# Inser the noise or correlation files not present in the event:
		if(length(out$noisefiles)==0){
			out$noisefiles=thisout$noisefiles
			}
		if(length(out$corfiles)==0){
			out$corfiles=thisout$corfiles
			}
		
		# Get the files closest to the specified time point, else pick out the first file of each type:
		if(length(out$noisefiles)>1){
			out$noisefiles = out$noisefiles[if(length(utim)>0) which.min(abs(out$noise_utim-utim)) else 1]
			}
		if(length(out$corfiles)>1){
			out$corfiles = out$corfiles[if(length(utim)>0) which.min(abs(out$cor_utim-utim)) else 1]
			}
		}
	
		
	########## Output ##########
	# Merge the noise files:
	out=unlist(out[c("noisefiles","corfiles")])
	##################################################
	##################################################
	}
