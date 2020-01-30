#*********************************************
#*********************************************
#' Create theoretical segmentation data for echosounder.
#'
#' @param event	The path to the event holding TSD files.
#' @param t		The time steps to generate theoretical segmentation for.
#'
#' @return
#'
#' @importFrom ncdf4 ncdim_def ncvar_def nc_create ncvar_def ncvar_put nc_close
#'
#' @export
#'
segmaskEchosounder <- function(event, t = "all") {
	
	if(identical(t, "all")) {
		# Get the number of pings:
		numt <- read.event(event = event, var = "numt")$numt
		t <- seq_len(numt)
	}
	
	if(is.list(event)) {
		event <- event$path
	}
	
	u <- UNIX_time(event)
	pingsIndt <- u$i000[endsWith(names(u$i000), "pings")]
	
	segFiles <- mapply(
		segmaskEchosounderOnePingsFile, 
		pingsIndt, 
		filePath = names(pingsIndt), 
		event = event
	)
	
	# Define the dimensions:
	dims <- list(
		segM = list(
			Name = c("sample", "beam", "ping"), 
			Unit = c("SampleIndex", "BeamIndex", "PingIndex")
		)
	)
	
	seg2NetCDF4(event)
}

segmaskEchosounderOnePingsFile <- function(t, filePath, event) {
	
	# Read the beams configuration data and extract the length and number of beams (assume simulated data with only one beam configuration):
	beams <- read.event(event, var = "beams", t = t, onestep = 2)
	lenb <- beams$lenb
	numb <- beams$numb
	
	utim <- read.event(event, var = "utim", t = t)
	
	# Read the school dynamics:
	schools <- read.event(event, var = labl.TSD("cs"), t = "all", drop = FALSE)
	schoolNR <- as.numeric(getSchoolIDs(getSchoolDirs(event)))
	dim(schoolNR) <- c(1, length(schoolNR))
	schools$schoolNR <- schoolNR
	
	# Get the segmentation masks of each ping of the file:
	segM <- lapply(t, segmaskEchosounderOnePing, event = event, schools = schools, beams = beams)
	# Unlist and set dimensions equal to the vbsc:
	segM <- unlist(segM)
	dim(segM) <- c(max(lenb), max(numb), length(t))
	
	# Create a list of TSD data:
	segData <- list(
		indt = t, 
		utim = utim$utim, 
		lenb = lenb, 
		numb = numb, 
		segM = segM
	)
	
	# Define the path to the netCDF4 file:
	con <- gsub(".pings$", ".seg", filePath)
	# Write the list to a -seg file:
	write.TSD(segData, con)
	
	# Return the file path:
	con
}

segmaskEchosounderOnePing <- function(t, event, schools, beams) {
	
	# Read the data of the ping:
	voxels <- read.event(event, var = "voxels", t = t)
	
	# Generate a segmentatnio mask for each school:
	masks <- lapply(seq_along(schools$schoolNR), segmaskEchosounderOnePingOneSchool, voxels = voxels, schools = schools, beams = beams)
	
	# Merge the masks, keeping only the largest integer value, assuming that the mixed categories are sorted later than the clean categories:
	mask <- do.call(pmax, masks)
	
	mask
}

segmaskEchosounderOnePingOneSchool <- function(schoolNR, voxels, schools, beams) {
	
	
	school <- lapply(schools, "[", , schoolNR)
	
	
	# Get the vector from each voxel to each school inside a given radius around the x,y of the first voxel:
	relativeX <- outer(voxels$psxx, school$psxS, "-")
	relativeY <- outer(voxels$psyx, school$psyS, "-")
	dist <- sqrt(relativeX^2 + relativeY^2)
	absoluteAngle <- atan2(relativeY, relativeX)
	
	# Get and subtract the voxel radius:
	#voxelRange <- (row(voxels$psxx) - 1) * beams$rres
	voxelRange <- array(seq_len(nrow(relativeX)) - 1, dim = dim(relativeX)) * beams$rres[1]
	voxelRadius <- voxelRange * tan(beams$bwtx[1]/2)
	distSubtracted <- dist - voxelRadius
	
	# Get and subtract the rotation angle:
	#oazS <- matrix(school$oazS, ncol = ncol(relativeX), nrow = nrow(relativeX), byrow = TRUE)
	oazS <- outer(array(0, dim = dim(voxels$psxx)), school$oazS, "+")
	relativeAngle <- absoluteAngle - oazS
	
	# Convert back to the x and y of the vector from a school to the voxel in the coordinate system of the school:
	distX <- distSubtracted * sin(relativeAngle)
	distY <- distSubtracted * cos(relativeAngle)
	
	# Get the vertical distance between the voxels and the school center:
	relativeZ <- outer(voxels$pszx, school$pszS, "-")
	# Get the radius of the ellipsoidal school at the vertical position of each voxel:
	#szxS <- matrix(school$szxS, ncol = ncol(relativeX), nrow = nrow(relativeX), byrow = TRUE)
	#szyS <- matrix(school$szyS, ncol = ncol(relativeX), nrow = nrow(relativeX), byrow = TRUE)
	#szzS <- matrix(school$szzS, ncol = ncol(relativeX), nrow = nrow(relativeX), byrow = TRUE)
	semiaxisX <- outer(array(0, dim = dim(voxels$psxx)), school$szxS / 2, "+")
	semiaxisY <- outer(array(0, dim = dim(voxels$psxx)), school$szyS / 2, "+")
	semiaxisZ <- outer(array(0, dim = dim(voxels$psxx)), school$szzS / 2, "+")
	radiusInVoxelDirection <- firstEllipseAxisFromSecond(
		y = relativeZ, 
		a = semiaxisX, 
		b = semiaxisZ
	)
	# Get factor to multiply the school by:
	schoolScale <- radiusInVoxelDirection / semiaxisX
	semiaxisXEffective <- schoolScale * semiaxisX
	semiaxisYEffective <- schoolScale * semiaxisY
	
	
	# Test for voxels inside the ellipsoids:
	mask <- distX^2 / semiaxisXEffective^2 + distY^2 / semiaxisYEffective^2 <= 1
	# Collapse across schools (of the same species / schoolNR):
	mask <- apply(mask, 1:2, max, na.rm = TRUE)
	# Ensure that the mask has non-negative values (-Inf generated if there are only NAs):
	mask[mask < 0] <- 0
	
	mask <- mask * schoolNR
}

firstEllipseAxisFromSecond <- function(y, a, b) {
	suppressWarnings(out <- sqrt(a^2 * (1 - y^2/b^2)))
	#out[is.na(out)] <- 0
	out
}

getSchoolDirs <- function(event) {
	files <- list.files(event, full.names = TRUE, pattern = "^school|School")
	dirs <- files[ file.info(files)$isdir ]
	dirs
}

getSchoolIDs <- function(schooldirs) {
	gsub("school|School", "", basename(schooldirs))
}
