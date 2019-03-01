#*********************************************
#*********************************************
#' Returns specifications of the acoustic instrument 'esnm' regarding types of beams. The function treats single time steps, except for the output "region".	
#'
#' @param data			A list of elements named according to the TSD file format holding the following information: "esnm" (may be given separately), "lenb", "", "", "", "", "", "", "", "", "", "", "".
#' @param var			The variables for which to get the specifications. "rect": get rectangular and non-rectangular beams. "offset": the angle offset and whether there is time offset. "region": The region over which to draw positions. "valid", "ind": The logical and indices of valid beams given \code{beamstypes}.
#' @param esnm			The name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param ind			The input indices of beams to use, which are updated by this funciton.
#' @param add			The angle in degrees to add on either side of the sonar volume.
#' @param beamstypes	The beam type to use, where rectangluar is 1 and non-rectangular is 2.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD strff zeros
#' @importFrom utils tail head
#'
#' @export
#' @rdname get.specs.esnm
#'
get.specs.esnm <- function(data, var=c("rect", "offset", "region", "valid", "ind"), esnm=NULL, ind=NULL, add=30, beamstypes=1){
			
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-09-06 - Clean version.
	
	##### Preparation and execution #####
	# Warning if multiple time steps are given, and not only "region" is requested:
	if(length(dim(data$dire)) > 1 && strff("off",var)){
		warning("Multiple time steps are given to the function 'get.specs.esnm', which only treats the first time steps when 'offset' is requested.")
	}
	
	# The name of the acoustical system may be given in 'data':
	if(length(esnm) == 0){
		esnm <- data$esnm[1]
	}
	

	# (1) Define which beams are considered circular and which are considrered rectangular, and set the dimensions of the output in both the rectangular and circular case:
	if(any(strff(c("val", "ind", "rec"),var))){
		# If not present in 'data', try to extract the number of beams
		if(length(data$numb)==0){
			numb <- sapply(data[c("dira", "dire", "lenb", "freq")],length)
			data$numb <- mean(numb[numb>0])
		}
		if(sonR_implemented(esnm, c("MBE"))){
			# Define the rectangular beams:
			nsplitbeams <- 2
			rect <- c(rep(TRUE, max(data$numb)-nsplitbeams), rep(FALSE, nsplitbeams))
		}
		else if(sonR_implemented(esnm, c("SBE"))){
			# Define the rectangular beams:
			rect <- rep(FALSE, length(data$freq))
		}
		else if(sonR_implemented(esnm, c("MBS","OFS"))){
			# Define the rectangular beams, and support multiple time steps but return only for one time step, thus using NROW which accepts both a vector and a matrix:
			rect <- rep(TRUE, NROW(data$freq))
	}
		else{
			# Define the rectangular beams:
			rect <- rep(TRUE, length(data$freq))
		}
	}
	else{
		rect <- NA
	}
				
	
	if(strff("off",var)){
		# (2) Wheter vessel time offset should be applied:
		hasOffset <- c("ms70")
		hasOffset <- strff(hasOffset,esnm)
			
		# (3) Angular offset to be applied when calculating the voxel edges. Note that the output from voxel.edge.TSD() is rotated by these angles and does not represent the true positions of the voxels. Therefore, all results of calculations based on the output from voxel.edge.TSD() must be rotated back to the original orientation using the negative of the engles in 'offset':
		if(sonR_implemented(esnm, c("MBE"))){
			offset <- list(by="x",ang=-pi/2)
		}
		else if(sonR_implemented(esnm, c("SBE"))){
			if(length(ind)>1 && length(ind[[2]])==1){
				beam <- ind[[2]]
			}
			else{
				beam <- 1
			}
			offset=list(by=c("z", "x"), ang=cbind(-data$dira[beam],-data$dire[beam]))
		}
		else if(sonR_implemented(esnm, c("OFS"))){
			if(length(data$bmmd)==0){
				vertical <- logical(NCOL(data$dire))
			}
			else{
				vertical <- head(data$bmmd, 1) == 2
			}
			ang <- zeros(length(vertical), 2)
			#ang[vertical,] <- cbind(head(data$dira, 1)[vertical], -pi/2)
			ang[vertical,] <- cbind(head(data$dira, 1)[vertical], pi/2 - head(data$dira, 1)[vertical])
			offset <- list(by=c("z", "x"), ang=ang)
		}
		else{
			offset <- list()
		}
	}
	else{
		hasOffset <- NA
		offset <- NA
	}
		
	# (4) Get the region over which to draw positions:
	if(strff("reg",var)){
		# 'add' is given in degrees:
		add <- add * pi/180
		if(sonR_implemented(esnm, c("MBS"))){
			region <- cbind(
				theta = c(min(data$dira) - add, max(data$dira) + add), 
				phi = c(min(data$dire) - add, max(data$dire) + add)
			)
		}
		else if(sonR_implemented(esnm, c("MBE"))){
			region <- cbind(theta=c(0,2*pi), phi=c(min(data$dire) - add,pi))
		}
		else if(sonR_implemented(esnm, c("SBE"))){
			region <- cbind(theta=c(0,2*pi), phi=c(pi - add,pi))
		}
		else if(sonR_implemented(esnm, c("OFS"))){
			if(length(data$bmmd)==0 || data$bmmd[1]==0){
				region <- cbind(theta=c(0,2*pi), phi=c(min(data$dire) - add,max(data$dire)+add))
			}
			else{
				region <- cbind(theta=c(0,2*pi), phi=c(pi/2-add,pi))
			}
		}
		else if(length(esnm)>0){
			stop(paste0("Acoustical instrument \"", esnm, "\", specified by 'esnm', not implemented"))
		}
		else{
			stop("No acoustical instrument specified by 'esnm'")
		}
	}
	else{
		region <- NA
	}
		
	# (5) Return the valid beams based on tha value of 'beamstypes':
	if(any(strff(c("val", "ind"), var))){
		valid <- list(which(rect), which(!rect))
		valid <- valid[sapply(valid,length)>0]
		valid <- unlist(valid[beamstypes])
		if(length(ind)>1){
			ind[[2]] <- intersect(ind[[2]], valid)
		}
		else{
			ind[[2]] <- valid
		}
	}
	else{
		ind=NA
	}
		
			
	##### Output #####
	# Add to the data:
	data$rect <- rect
	data$ind <- ind
	data$hasOffset <- hasOffset
	data$offset <- offset
	data$region <- region
	data$valid <- valid
	data
	#output=list(esnm=esnm, rect=rect, ind=ind, hasOffset=hasOffset, offset=offset, region=region, valid=valid)
}
