#*********************************************
#*********************************************
#' Transforms from 'sgsc' to positions, in a similar way as used in pplot3d.TSD(). Used in read.event() for clustering 'sgsc'.
#'
#' @param filelist  is a vector of all file names.
#' @param filesind  is a vector of the indexes of the segmentation files in the list of files.
#' @param segpar  is a list of elements named "bwGp", "lsth"/"rlst", "usth"/"rust", or "sgth"/"sgt0" specifying the parameters of the segmentation data to read.
#' @param TIME  is the list returned from UNIX_time().
#' @param tlist  is the list of time indexes to be read for each file.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSDs
#'
#' @export
#' @rdname read.event_read_sgsc
#'
read.event_read_sgsc<-function(filelist, var, filesind, segpar, TIME, tlist, merge, msg=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-08-28 - Clean version.
	# Update: 2012-08-29 - Changed to select the value of the segmentation parameters that is closest to the one given in 'segpar'.
	# Update: 2012-10-02 - Changed according to the change in sgPM, ensuring that this matrix is coordinated with the sorted list of files.
	# Update: 2012-10-10 - Changed according to the change in sgPM.
	# Update: 2012-10-14 - Fixed bug when inputing segpar as a numeric index.
	# Last: 2013-10-09 - Included removeDuplicated_tlist() in read.event_read_sgsc() in order to avoid warnings associated with duplicated time steps before the selection of segmentation files is done.
	########### DESCRIPTION: ###########
	# Transforms from 'sgsc' to positions, in a similar way as used in pplot3d.TSD(). Used in read.event() for clustering 'sgsc'.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---filelist--- is a vector of all file names.
	# ---filesind--- is a vector of the indexes of the segmentation files in the list of files.
	# ---segpar--- is a list of elements named "bwGp", "lsth"/"rlst", "usth"/"rust", or "sgth"/"sgt0" specifying the parameters of the segmentation data to read.
	# ---TIME--- is the list returned from UNIX_time().
	# ---tlist--- is the list of time indexes to be read for each file.
	

	##################################################
	##################################################
	if(!is.numeric(segpar)){
		# 'segpar' may be given with more user friendly variables:
		if("h" %in% names(segpar) && !"bwGp" %in% names(segpar)){
			segpar$bwGp=segpar$h
			}
		if("beta0" %in% names(segpar) && !"lsth" %in% names(segpar)){
			segpar$lsth=segpar$beta0
			}
		if("beta1" %in% names(segpar) && !"usth" %in% names(segpar)){
			segpar$usth=segpar$beta1
			}
		if("rbeta0" %in% names(segpar) && !"rlst" %in% names(segpar)){
			segpar$rlst=segpar$rbeta0
			}
		if("rbeta1" %in% names(segpar) && !"rust" %in% names(segpar)){
			segpar$rust=segpar$rbeta1
			}
		if(any(c("c","alpha") %in% names(segpar)) && !"sgth" %in% names(segpar)){
			segpar$sgth=segpar$c
			}
		if("c0" %in% names(segpar) && !"sgt0" %in% names(segpar)){
			segpar$sgt0=segpar$c0
			}
		
		# Read the parameters of the segmentation, stored in the file names of the segmentation files:
		thissegpar=sgPM(filelist[filesind])
		if(length(thissegpar)==0){
			return()
			}
		thissegpar=cbind(thissegpar,SFNR=seq_len(nrow(thissegpar)))
		
		# Remove the variables that are not present in the files:
		segpar=segpar[intersect(names(segpar),colnames(thissegpar))]
		
		# Select the segmentation file(s) based on the information given in 'segpar':
		segfile=seq_along(filelist[filesind])
		
		
		############ Pick out the closest value, and not requiring to match exactly: ##############
		if(length(segpar)>0){
			if(is.list(segpar)){
				warn=NULL
				segfile=!logical(nrow(thissegpar))
				
				for(i in seq_along(segpar)){
					if(is.na(segpar[[i]])){
						segfile = segfile & is.na(thissegpar[,names(segpar[i])])
						}
					else{
						difference=abs(segpar[[i]]-thissegpar[,names(segpar[i])])
						minvalue=max(sqrt(.Machine$double.eps),min(abs(thissegpar[,names(segpar[i])]),na.rm=TRUE))*sqrt(.Machine$double.eps)
						if(any(difference < minvalue,na.rm=TRUE)){
							segfile = segfile & difference < minvalue
							}
						else if(!all(is.na(thissegpar[,names(segpar[i])]))){
							warn=c(warn,paste0("The closest value of ", names(segpar[i]), " chosen (", signif(thissegpar[,names(segpar[i])][which.min(difference)], 5), " vs input value ", segpar[i], ")"))
							segfile = segfile & difference < min(difference,na.rm=TRUE) + minvalue
							}
						}
					}
				# Issue the warnings if any:
				if(length(warn)>0){
					warning(paste(warn,collapse="\n",sep=""))
					}
				# Convert to numeric indices:
				segfile=which(segfile)
				
				# Check if the selection of segmentation files is unique:
				if(nrow(unique(thissegpar[segfile,,drop=FALSE]))>1){
					##### THE FUNCITON DOES NOT ACCEPT PINGS DISTRIBUTED OVER MULTIPLE FILES. THIS NEEDS FIXING!!!!!!!!! #####
					warning(paste("Non-unique selection of segmentation data using the list 'segpar' of segmentation parameters. The following segmentation file nrs were selected:",paste(segfile,collapse=" "),sep="\n"))
					}
				# Use the file indexes stored in the sgPM-file:
				segfile=thissegpar[segfile,"ORDR"]
				}
			else if(is.numeric(segpar)){
				segfile=intersect(segpar,seq_along(filelist[filesind]))
				}	
			else{
				warning("'segpar' must be a list of elements named \"bwGp\", \"lsth\", \"usth\" (or \"rlst\", \"rust\"), \"sgth\", \"code\" or other parameters present in the segmentation files, specifying which of the segmentation files to read. Else all files are read and returned as a list.")
				}
			}
		else if(length(filelist[filesind])>1){
			warning("No selection of segmentation files applied using 'segpar'. All segmentation files were included.")
			}
		
		# Read the segmentation file(s):
		segfile=as.numeric(segfile)
		segparinfo = paste(c(rbind(names(segpar),": ", unlist(segpar),"\n")), collapse="")
		if(msg){
			cat(ngettext(length(segfile), "The following segmentation file selected", "The following segmentation files selected"),"\n",paste(seq_along(segfile),filelist[filesind][segfile],collapse="\n"),"\n", "for the following segmentation parameters:\n", segparinfo)
		}
		
		suppressWarnings(c(read.TSDs(filelist[filesind][segfile], t=tlist[filesind[segfile]], var=var, dimension=TRUE, merge=merge, clean=FALSE, indt=FALSE, msg=FALSE, addNvar=TRUE),list(sgPM=thissegpar, segfilenr=segfile, segfile=filelist[filesind][segfile])))
		}
	# Else the segmentation file numbers are already given (segpar numeric):
	else{
		segfile=segpar
		segparinfo = paste("Number of the file in the sorted list of files: ", segpar, "\n", sep="")
		if(msg){
			cat(ngettext(length(segfile), "The following segmentation file selected", "The following segmentation files selected"),"\n",paste(seq_along(segfile),filelist[filesind][segfile],collapse="\n"),"\n", "for the following segmentation parameter:\n", segparinfo)
		}
		
		suppressWarnings(c(read.TSDs(filelist[filesind][segfile], t=tlist[filesind[segfile]], var=var, dimension=TRUE, merge=merge, clean=FALSE, indt=FALSE, msg=FALSE, addNvar=TRUE),list(segfilenr=segfile, segfile=filelist[filesind][segfile])))
		}
	##################################################
	##################################################
	}
