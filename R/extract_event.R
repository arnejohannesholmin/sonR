#*********************************************
#*********************************************
#' Extracts events from a directory of raw files.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param t  is a two element vector of ftim values between which rawfiles are extracted from a directory (see the description of 'event').
#' @param rawevent  specifies a directory with raw files which are copied to 'event' and converted to TSD filess before running the segmentation.
#' @param ow  is FALSE to not overwrite existing data when copying data from 'rawevent'.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param toTSD  is false to disable converting to TSD files.
#' @param last  is FALSE to skip the last raw file, useful if this file is not complete.
#' @param temp  is TRUE to only return the temporary directory if raw files already exist and toTSD==TRUE, in which case the existing and temporary tsd files are not merged into the existing.
#' @param ... are parameters passed to EKRaw2TSD().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ftim2utim utim2ftim
#' @importFrom tools file_ext
#' @importFrom utils tail head
#'
#' @export
#' @rdname extract_event
#'
extract_event <- function(event, rawevent, t, ow=TRUE, dir.data=NULL, toTSD=TRUE, last=TRUE, temp=FALSE, ...){
	
	############### LOG: ###############
	# Start: 2015-10-05 - Clean version.
	
	# Repeat the t if only of length 1:
	if(length(t) == 1){
		t = rep(t,2)
		}
	# Allow for strings of the type "yyyy-mm-dd HH:MM:SS" or other ways of writing the time:
	t = utim2ftim(ftim2utim(t))
	# Get start times of the files:
	l = list.files(rawevent[1], full.names = TRUE)
	l = l[tolower(file_ext(l)) == "raw"]
	basel = basename(l)
	
	atD = sapply(gregexpr("D", basel, fixed = TRUE), head, 1)
	atT = sapply(gregexpr("-T", basel, fixed = TRUE), head, 1)
	D = substr(basel, atD+1, atD+1+7)
	T = substr(basel, atT+2, atT+2+5)
	ftim1 = as.numeric(paste0(D, T))
	utim1 = ftim2utim(ftim1)
	lastutim = tail(ftim1, 1) + if(length(ftim1)>1) diff(tail(ftim1,2)) else 0
	lastutim = min(unclass(Sys.time()), lastutim)
	utim2 = c(utim1[-1],lastutim)
	ftim2 = utim2ftim(utim2)
	
	# Pick out the files:
	valid = which(ftim2 >= as.numeric(t[1]) & ftim1 <= as.numeric(t[2]))
	# Add the file before and after the sequence:
	if(length(valid)>0){
		if(!last){
			valid = valid[length(valid)]
			}
		}
	l = l[valid]
	
	# Generate the event directory in which to put the raw files before converting to TSD:
	if(length(event)>0){
		if(length(event)==3){
			event = generate.event(cruise=event[1], event=event[2], esnm=event[3], dir.type="raw", dir.data=dir.data)
			}
		else{
			# Define the raw and tsd target directory
			event = file.path(dirname(event[1]), "raw")
			dir.create(event, recursive=TRUE)
			}
		}
	event_raw = file.path(dirname(event[1]), "raw")
	event_tsd = file.path(dirname(event[1]), "tsd")
	
	# If files already extist in the event_raw, create a temporary event with the new files:
	if(toTSD){
		# Get the existing raw files:
		rawfiles = list.files(event_raw, "^.*[:.:](raw)$")
		newrawfiles = !(basename(l) %in% rawfiles)
		
		tempevent_tsd = NULL
		# If any raw files already exist in the target directory:
		if(length(rawfiles)>0){
			tempevent_raw = file.path(dirname(event_raw), "temp_extract_event", "raw")
			tempevent_tsd = file.path(dirname(tempevent_raw), "tsd")
			dir.create(tempevent_raw, recursive=TRUE)
			
			# If there are new raw files to copy to the target directory:
			if(sum(newrawfiles)>0){
				# Copy the raw files to the tempevent_raw and convert to TSD files:
				
				cat("The following raw files will be copied to the temporary event \"",tempevent_raw,"\":\n",paste0(seq_along(l[newrawfiles]), ". ",l[newrawfiles], "\n", collapse=""))
				file.copy(l[newrawfiles], file.path(tempevent_raw, basename(l[newrawfiles])))
				# Convert to TSD files:
				EKRaw2TSD(tempevent_raw, dir.data="", ow=ow, ...)
				tempevent_tsd=file.path(dirname(tempevent_raw), "tsd")
				# Merge the raw and tsd events:
				file.copy(list.files(tempevent_raw, full.names=TRUE), event_raw)
				unlink(tempevent_raw, recursive=TRUE)
				if(!temp){
					combine.events(c(event_tsd, tempevent_tsd), event_tsd, cruise=NULL, esnm=NULL, dir.data=NULL, ctd=1)
					unlink(tempevent_tsd, recursive=TRUE)
					t = read.event(event_tsd, var="indt", t=t)$indt
					}
				else{
					t = read.event(tempevent_tsd, var="indt", t=t)$indt
					}
				}
			# If no new raw files will be converted:
			else{
				t = read.event(event_tsd, var="indt", t=t)$indt
				}
			}
		# If no raw files exist in the target directory:
		else if(length(rawfiles)==0){
			cat("The following raw files will be copied to the event \"",event_raw,"\":\n",paste0(seq_along(l), ". ",l, "\n", collapse=""))
			file.copy(l, file.path(event_raw, basename(l)))
			EKRaw2TSD(event_raw, dir.data="", ow=ow, ...)
			t = read.event(event_tsd, var="indt", t=t)$indt
			}
		list(event=event_tsd, t=t, tempevent=tempevent_tsd)
		}
	else{
		# Copy the raw files to the rawevent and convert to TSD files:
		cat("The following raw files will copied to the event \"",event_raw,"\":\n",paste0(seq_along(l), ". ",l, "\n", collapse=""))
		file.copy(l, event_raw)
		event_tsd
		}
	}
