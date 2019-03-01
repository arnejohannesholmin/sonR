#*********************************************
#*********************************************
#' (Internal) Merge segmentation events, used in SX90.segment.event() and extract_event().
#'
#' @param events  is a list of beam configuration data.
#' @param event  is a list of beam configuration data.
#' @param cruise  is a list of beam configuration data.
#' @param esnm  is a list of beam configuration data.
#' @param dir.data  is a list of beam configuration data.
#' @param ctd  is a list of beam configuration data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD combine.TSD read.TSD read.TSDs zeropad
#' @importFrom tools file_ext
#'
#' @export
#' 
combine.events <- function(events, event, cruise=NULL, esnm=NULL, dir.data=NULL, ctd=1){
	if(length(events)<2){
		stop("At least one event must be given")
		}
	event = generate.event(event=event, cruise=cruise, esnm=esnm, dir.type="tsd", dir.data=dir.data)
	# Get the paths to the files of type "pings", "beams", "ctd", "vessel", and "rawvessel":	
	fileslist = lapply(events, list.files, full.names=TRUE)
	files = unlist(fileslist)
	ext = file_ext(files)
	pingsfiles = files[ext=="pings"]
	vesselfiles = files[ext=="vessel"]
	ctdfiles = files[ext=="ctd"]
	rest = setdiff(files, c(pingsfiles, vesselfiles, ctdfiles))
	# Get the variables in the rest of the files:
	varsinrest = read.TSDs(rest, var="", t=1, clean=FALSE, merge=FALSE, header=TRUE)
	varsinrest = varsinrest[names(varsinrest)=="labl"]
	beamsfiles = rest[sapply(varsinrest, function(x) "freq" %in% x)]
	rawvesselfiles = rest[sapply(varsinrest, function(x) "irzv" %in% x)]
	
	# Copy the seg-files to the new event:
	get_sfnr = function(x, replace=NULL){
		atsfnr = unlist(gregexpr("sfnr_", x, fixed=TRUE)) + 5
		atseg = unlist(gregexpr(".seg", x, fixed=TRUE)) - 1
		if(length(replace)>0){
			paste0(substr(x, 1, atsfnr-1), replace, substring(x, atseg+1))
			}
		else{
			unlist(substr(x, atsfnr, atseg))
			}
		}
	
	segfileslist = lapply(fileslist, function(x) x[ext=="seg"])
	if(all(sapply(segfileslist, length)>0)){
		segfiles = unlist(segfileslist)
		sfnr = lapply(segfileslist, get_sfnr)
		# Get maximum sfnr in events[1], and insert sfnr preceding this maximum value in events[2]:
		max1 = max(as.numeric(sfnr[[1]]))
		sfnr[[2]] = zeropad(seq_along(sfnr[[2]])+max1, nchar(sfnr[[2]][1]))
		targetsegfiles = get_sfnr(segfiles, replace=unlist(sfnr))
		file.copy(segfiles, file.path(event, basename(targetsegfiles)))
		}
	
	# Copy the pings-files to the new event:
	file.copy(pingsfiles, event)
	
	# Merge beams files if these have more than one time step, otherwise copy only the first file to the new event:
	if(read.TSD(beamsfiles[1], var="numt")$numt>1){
		combine.TSD(beamsfiles, dir=event, indt=FALSE, reserve=FALSE)
		}
	else{
		file.copy(beamsfiles[1], file.path(event,basename(beamsfiles[1])))
		}
	
	# Merge vessel files:
	combine.TSD(vesselfiles, dir=event, indt=FALSE, reserve=FALSE)

	# Merge rawvessel files:
	combine.TSD(rawvesselfiles, dir=event, indt=FALSE, reserve=FALSE)

	# Copy the first CTD file to the new event:
	file.copy(ctdfiles[ctd], file.path(event,basename(ctdfiles[1])))
	}
