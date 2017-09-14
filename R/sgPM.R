#*********************************************
#*********************************************
#' Reads the ordered matrix of segmentation parameters from the files given in 'segfiles', or attempts to read the variable 'sgPM' directly if 'segfiles' is a directory. If 'segfiles' is a directory, and the variable 'sgPM' is not present, a file is written with the matrix of segmentation parameters ("sgPM.tsd").
#'
#' @param segfiles  is a string or vector of strings representing the path to the files or the directory from which to obtain the UNIX_time information, either by reading the UNIX_time file or by generating new information.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSD strff write.TSD
#' @importFrom tools file_ext

#'
#' @export
#' @rdname sgPM
#'
sgPM<-function(segfiles){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-09-19 - Clean version.
	# Update: 2012-10-02 - Changed to return 9 columns, where column 8 is the size of the file and column 9 is the file index in the sorted list of files. Column 9 is used to coordinate reading and the sgPM.
	# Update: 2012-10-10 - Changed to return 11 columns, adding the relative lower and upper schooling threshold.
	# Last: 2013-07-16 - Changed to the new format of variable names and values enclosed in <>. This can accept any variable name.
	########### DESCRIPTION: ###########
	# Reads the ordered matrix of segmentation parameters from the files given in 'segfiles', or attempts to read the variable 'sgPM' directly if 'segfiles' is a directory. If 'segfiles' is a directory, and the variable 'sgPM' is not present, a file is written with the matrix of segmentation parameters ("sgPM.tsd").
	########## DEPENDENCIES: ###########
	# read.TSD(), write.TSD()
	############ VARIABLES: ############
	# ---segfiles--- is a string or vector of strings representing the path to the files or the directory from which to obtain the UNIX_time information, either by reading the UNIX_time file or by generating new information.
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Vector of the variable names formerly used in the file names of segmentation files:
	used=c("bwG1", "bwG2", "bwG3", "bwGp", "lsth", "usth", "rlst", "rust", "sgth", "sgt0", "misM", "MSVS", "code", "mdf1", "mdf2", "mdf3", "mdft", "nsTh", "dBan", "dBb1", "dBb2", "dBb3", "dBbs", "smty", "sfnr")
	
	# Get the directory of the segmentation files:
	if(length(segfiles)==1 && isTRUE(file.info(segfiles)$isdir)){
		dir=segfiles
		}
	else{
		dir=dirname(segfiles[1])
		}
	# Extract the segmentation files if 'segfiles' is a directory:
	if(length(segfiles)<2 && !strff("seg", file_ext(segfiles))){
		segfiles=list.files_caseInsensitive(dir,full.names=TRUE,recursive=TRUE)
		ext = file_ext(segfiles)
		segfiles=segfiles[ext=="seg"]
		}
	
	# Return an empty list if there are no valid segfiles:
	if(length(segfiles)==0){
		return(list())
		}
	
	# Sort the files:
	segfiles=sort(segfiles)
	
	# Identify the tsd directory, if present:
	tsddir=gregexpr("/tsd",dir)[[1]]
	if(tsddir[1]>0){
		attsd=max(tsddir)
		dir=substr(dir,1,attsd+3)
		}
	
		
	##### Execution and output #####
	# Declare the output:
	out=NULL
	
	# Look for the sgPM.tsd file:
	sgPMfile=file.path(dir,"sgPM.tsd")
	# If the file "sgPM.tsd" exists, extract the information:
	if(file.exists(sgPMfile)){
		# Read the "sgPM.tsd" file:
		sgPM=read.TSD(sgPMfile,var=c("sgPM","sgPN"),dimension=TRUE,header=FALSE,t="all")
		out=sgPM$sgPM
		if(NROW(out)!=length(segfiles)){
			out=NULL
			}
		else{
			colnames(out)=sgPM$sgPN
			}
		}
	
	if(length(out)==0 && length(segfiles)>0){
		# If 'segfiles' is given as a vector of file names, check that the number of files equals the number of rows of 'out':
		if(length(segfiles)!=NROW(out)){
			# Read the parameters of the segmentation, stored in the segmentation files:
			preout=list()
			for(i in seq_along(segfiles)){
				
				thesevar=read.TSD(segfiles[i],var=used,header=FALSE,info=FALSE)
				if(length(thesevar$bwGp)==3){
					thesevar$bwGp_x=thesevar$bwGp[1]
					thesevar$bwGp_y=thesevar$bwGp[2]
					thesevar$bwGp_z=thesevar$bwGp[3]
					thesevar$bwGp=thesevar$bwGp[1]
					}
				if(length(thesevar$mdft)==3){
					thesevar$mdf1=thesevar$mdft[1]
					thesevar$mdf2=thesevar$mdft[2]
					thesevar$mdf3=thesevar$mdft[3]
					thesevar$mdft=thesevar$mdft[1]
					}
				
				preout[[i]]=thesevar
				# Get the size of the file:
				preout[[i]]$SIZE=file.info(segfiles[i])$size
				}
			# Get all the variable names contained in the file names:	
			allvar=unique(unlist(lapply(preout,names)))
			
			# Convert to a data frame of the files along the first dimension and the variable along the second dimension (we use data frame here to allow for different types of the variables in the future (if the read.TSD() and write.TSD() are expanded to support writing data frames) and to ease the reordering of the columns below. However, data frames are slow and should be generally avoided.):
			out=matrix(NA,nrow=length(segfiles),ncol=length(allvar)+1)
			colnames(out)=c(allvar,"ORDR")
			# Insert the data in the output matrix:
			for(i in seq_along(preout)){
				namesi=names(preout[[i]])
				for(j in seq_along(namesi)){
					out[i,namesi[j]]=preout[[i]][[j]][1]
					}				
				}
			
			# Add the indexes of the sorted files:	
			out[,"ORDR"]=seq_along(segfiles)
			allvar=c(allvar,"ORDR")
			
			# Reorder the columns of 'out':
			present=intersect(allvar,used)
			new=setdiff(allvar,used)
			out=out[,c(present,new),drop=FALSE]
			
			# Order and transform back to matrix:
			ORDR=do.call(order,as.data.frame(out[,present,drop=FALSE]))
			out=out[ORDR,,drop=FALSE]
			
			# Write the matrix to file:
			#cat("Writing segmentation parameter file\n")
			write.TSD(list(sgPM=out,sgPN=colnames(out)),header=list(dtyp=c("doub","char")),sgPMfile,dimension=TRUE)
			}
		}
	
	out
	##################################################
	##################################################
	}
