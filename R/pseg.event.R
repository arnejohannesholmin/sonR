#*********************************************
#*********************************************
#' Segment the original fish positions as the voxels that contain at least the fraction 'alpha' of the expected number of fish, calculated as the volume of the voxels times the packing density:
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param t  is either the number of the ping to be treated, as listed from 1 to the number of pings in the event, or the time point given as a string "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF". Only one time step alowed!
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event.
#' @param alpha  is the significance level of the hypothesis testing of H0: school not present in the voxel, against H1: school present in the voxel. To return non-thresholded data, use alpha=NULL.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param var  is "vbsc" for segmenting volume backscattering data, and "posf" for segmenting original fish position data.
#' @param dens  is the packing density of the school, used when var="posf". These values are attemptedly read in the event, if not given as input to the function.
#' @param segdir  is an optional string giving the name of the folder in which to put the segmentation files. The folder will be put in the directory of the event. segdir==NULL results in segdir="seg0".
#' @param adds  is an optional list of variables overriding the variables in read from the event.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw soundbeam_range
#' @importFrom TSD arr.ind2ind read.TSD strff write.TSD zeros
#' @importFrom utils tail head
#'
#' @export
#' @rdname echoIBM.segment.event
#'
pseg.event <- function(event=1, t=1, cruise=2009116, alpha=NULL, filesize=3e8, dens=NULL, segdir=NULL, adds=list(), TOV=0, startn=1, ...){
	
	########## Preparation ##########
	# Get the directory of the event:
	event <- event.path(event=event, cruise=cruise)$event
	# Get the time indices and the formated time strings corresponding to 't':
	time <- read.event(event=event, var=c("indt", "time"), t=t)
	
	# The number of time steps:
	numt <- length(time$indt)
	# Read the beam configuration:
	beams <- read.event(event=event, var="beams")
	
	# 'totallength' limits the size of the files:
	totallength <- zeros(length(alpha))
	# Define the header used when writing the data:
	newfile <- !logical(length(alpha))
	
	# Create the directory to put the segmentation files into:
	if(length(segdir)==0){
		dir.out <- file.path(event, c("seg", "seg0")[startsWith("p", tolower(var[1])) + 1])
		}
	else{
		dir.out <- file.path(event, segdir)
		}
	if(!file.exists(dir.out)){
		suppressWarnings(dir.create(dir.out))
		}
	
	add <- ""
	
	
	########## Execution and output ##########
	# Segment the original fish positions as the voxels that contain at least the fraction 'alpha' of the expected number of fish, calculated as the volume of the voxels times the packing density:
		
	# Define the name of the segmentation file:
	header <- list(dtyp=c("floa", "doub", "floa", "floa", "floa", "floa", "doub", "floa"))
	
	# Read the packing density if not given as input to the function:
	if(length(dens)==0){
		dens <- read.event(event=event, t=time$indt, var="scpd")$scpd
		if(length(dens)==0){
			stop("The packing density must be given or present in a TSD file as a variable named \"scpd\"")
			}
		}
	# Get the radial, azimuth, and elevation edges in the spherical coordinate system with the sonar at the origin:
	J <- max(beams$lenb)
	I2 <- length(unique(beams$freq))
	I1 <- length(beams$freq)/I2
	numb <- length(beams$lenb)
	# Define the dimensions of the conical voxels (in the spherical coordinate system), and use these to define the vectors used in findInterval():
	r <- soundbeam_range(beams)
	
	utheta <- unique(beams$dira)
	upptheta <- TRUE
	if(utheta[1]>utheta[2]){
		upptheta <- FALSE
		}
	utheta <- sort(utheta)
	dtheta <- diff(utheta)[1]
	theta <- c(utheta[1]-dtheta/2, utheta+dtheta/2)
	
	uphi <- unique(beams$dire)
	uppphi <- TRUE
	if(uphi[1]>uphi[2]){
		uppphi <- FALSE
		}
	uphi <- sort(uphi)
	
	dphi <- diff(uphi)[1]
	phi <- c(uphi[1]-dphi/2, uphi+dphi/2)
	# Get the voxel volumes:
	volx <- volx.TSD(beams)
					
	# Read the vessel dynamics:
	vessel <- read.event(event=event, var="vessel", t=time$indt, TOV=TOV)
		
	cat("Segmenting original fish positions for time steps:\n")
	# Run through the time steps:
	for(i in seq_len(numt)){
		cat(time$indt[i], " of ", numt, "\n", sep="")
		# Read the fish position data:
		data <- read.event(event=event, var=c("psxf", "psyf", "pszf", "indt", "utim"), t=time$indt[i], cs="v")
		
		# Rotate into the coordinate system of the sonar:
		pos_sph <- rotate3D(cbind(data$psxf-vessel$psxv[i], data$psyf-vessel$psyv[i], data$pszf-vessel$pszv[i]-beams$psze), by="z", ang=vessel$rtzv[i], sph.out=TRUE)
		# Identify which voxels the fish are located in:
		fishr <- findInterval(pos_sph[, 1], r)
		if(upptheta){
			fishtheta <- findInterval(pos_sph[, 2] %% (2*pi), theta)
			}
		else{
			fishtheta <- length(theta) - findInterval(pos_sph[, 2] %% (2*pi), theta)
			}
		if(uppphi){
			fishphi <- findInterval(pos_sph[, 3], phi)
			}
		else{
			fishphi <- length(phi)-findInterval(pos_sph[, 3], phi)
			}
		vox <- TSD::arr.ind2ind(cbind(fishr, fishtheta, fishphi), c(J, I1, I2))
		
		# Calculate the expected number of fish in ecah voxel
		expected <- dens[i]*volx
		# Get the number of fish in each voxel:
		nfish <- table(vox)
		nfish <- cbind(as.numeric(names(nfish)), nfish)
		
		for(thisa in seq_along(alpha)){
			thissgs0 <- as.integer(nfish[nfish[, 2]>alpha[thisa]*expected[nfish[, 1]], 1])
			if(newfile[thisa]){
				segfile <- echoIBM.get.segfilename(n=length(thisa), event=dir.out[1], add=paste(add, time$indt[i], sep="_"), startn=startn)
				sfnr <- segfile[[2]]
				segfile <- segfile[[1]]
				# Write the data to file, adding the parameters of the segmentation, the bandwidth of the Gaussian kernel smoothing of the probability "bwGp" ('h'), the lower schooling threshold value "lsth" ('beta0'), the upper schooling threshold value of the prior distribution "usth" ('beta1'), and the segmentation threshold value "sgth" ('alpha'):
				suppressWarnings(bytes <- TSD::write.TSD(con=segfile[[1]], header=header, x=c(list(indt=time$indt[i], mtim=time$mtim[i]), beams[c("numb", "lenb", "freq")], list(sgs0=thissgs0, Ttvl=sum(volx[thissgs0]), sgt0=alpha[thisa]), segfile[2]), numt=1, reserve=numt-1, keep.null=TRUE, ...))
				}
			else{
				# Write (append) the data to file:
				suppressWarnings(bytes <- TSD::write.TSD(con=segfile, header=header, x=list(indt=time$indt[i], mtim=time$mtim[i], numb=NULL, lenb=NULL, freq=NULL, sgs0=thissgs0, Ttvl=sum(volx[thissgs0]), sgt0=NULL), numt=1, append=TRUE, keep.null=TRUE, ...))
				
				}
			totallength[thisa] <- totallength[thisa]+bytes
			if(totallength[thisa]>filesize){
				newfile[thisa] <- TRUE
				totallength[thisa] <- 0
				}
			else{
				newfile[thisa] <- FALSE
				}
			}
		}
	
	# Update the UNIX_time file, since time steps has been appended to the segmentation file(s):
	TSD::UNIX_time(event=event, file=TRUE, fresh=TRUE)
	
	invisible(totallength)
	##################################################
	##################################################
	}
