#*********************************************
#*********************************************
#' Read per-ping-output file from Profos (one line per ping per school).
#'
#' @param x					The Profos per ping file.
#' @param currentSpeed		A vector of the speed of the current at the observation, one value per ping or one value for all pings.
#' @param currentAngle		A vector of the angle of the current at the observation, in degrees clockwise from North, one value per ping or one value for all pings.
#' @param compensateCurrent	Logical: If TRUE (the default), compensate for the current.
#' @param sparfact			A factor to multiply the \code{spar} used in \code{\link{smooth.spline}} when smoothing the school positions.
#' @param spar0				The value of \code{spar} used if application of the \code{sparfact} results in negative values for \code{spar}.
#' @param ...				Arguments passed to \code{\link{lowess}}.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SoDA geoXY
#' @importFrom TSD ang2rot papply
#'
#' @export
#' @rdname readProfosPP
#'
readProfosPP <- function(x, currentSpeed=NULL, currentAngle=NULL, compensateCurrent=TRUE, sparfact=1, spar0=0.5, byId=TRUE, rm.duplicates=FALSE, filter=NULL, cores=1, nrows=-1, smooth=TRUE, ...){
	
	# TODO:
	#- Plot degrees 
	#- Add vertical lines for start of the net
	#
	#
	#Before/After
	#- Mean sv
	#- Polarisation
	#- Depth
	#- Speed
	
	# Function for subtracting the current if given:
	subtractCurrent <- function(x, vars){
		# Compensate for current, by subtracting the current speed multiplied by the time from the first time in each cartesian direction:
		if(!any(is.na(x$currentSpeed))){
			# Convert to m/s:
			currentSpeed <- currentSpeed * 1852 / 3600
			# Convert to radians on the unit circle:
			currentAngle <- currentAngle * pi / 180
			currentAngle <- pi/2 - currentAngle
			# Get speed in x and y:
			currentSpeed_x <- cos(currentAngle) * currentSpeed
			currentSpeed_y <- sin(currentAngle) * currentSpeed

			deltaSec <- x$UNIXtime - x$UNIXtime[1]
			x[[vars[1]]] <- x[[vars[1]]] - currentSpeed_y * deltaSec
			x[[vars[2]]] <- x[[vars[2]]] - currentSpeed_x * deltaSec
		}
		
		x
	}
	
	# Function for extracting diff, heading and speed of consecutive time steps, for smoothed or unsmoothed data (through the 'type' parameter):
	getDiffHeadingSpeed <- function(x, var="Center", type=""){
		
		diffWithNA <- function(x){
			c(diff(x), NA)
		}
		getSpeed <- function(dx, dy, dt){
			sqrt(dx^2 + dy^2) / dt
		}
		
		# Treat variable names:
		if(nchar(type)>0){
			type <- paste0("_", type)
		}
		vars <- paste0(var, ".", c("x", "y"), type)
		varsDiff <- paste0(var, ".", "d", c("x", "y"), type)
		varHeading <- paste0(var, ".", "heading", type)
		varSpeed <- paste0(var, ".", "speed", type)
		
		# Compensate for the current:
		x_comp <- subtractCurrent(x, vars)
			
		# Get diffs in x, y and time:
		d <- setNames(lapply(x_comp[c(vars, "UNIXtime")], diffWithNA), varsDiff)
		# Get heading from atan2 of y and x:
		h <- setNames(list(do.call(atan2, unname(d[2:1]))), varHeading)
		# Get speed from Eucledian distance of x and y, divided by diff of time:
		s <- setNames(list(do.call(getSpeed, unname(d))), varSpeed)
		
		cbind(x, d[1:2], h, s)
	}
	
	# Smooth one stretch of data (given by the column Id):
	smooth.spline_xy <- function(x, sparfact=NULL, spar0=0.5, var=c("Center.x", "Center.y"), xynames=c("Center.x_spline", "Center.y_spline")){
		# Smooth first to get the computed spar:
		fit_x <- smooth.spline(x[[var[1]]] ~ x$UNIXtime)
		fit_y <- smooth.spline(x[[var[2]]] ~ x$UNIXtime)
		#plot(x$Center.x, x$Center.y, type="o")
		#points(fit_x$y, fit_y$y, type="o", col=2)
		#print(c(fit_x$spar, fit_y$spar))
		if(length(sparfact)){
			# Multiply the resulting spar with sparfact, or if negative use the spar0:
			if(any(fit_x$spar < 0, fit_y$spar < 0)){
				sparx <- spar0
				spary <- spar0
			}
			else{
				sparx <- sparfact * fit_x$spar
				spary <- sparfact * fit_y$spar
			}

			fit_x <- smooth.spline(x[[var[1]]] ~ x$UNIXtime, spar=sparx)
			fit_y <- smooth.spline(x[[var[2]]] ~ x$UNIXtime, spar=spary)
		}
		#plot(x$Center.x, x$Center.y, type="o")
		#points(fit_x$y, fit_y$y, type="o", col=2)
		y <- setNames(list(fit_x$y, fit_y$y), xynames)
		x <- cbind(x, y)
		return(x)
	}
	
	# Smooth one stretch of data (given by the column Id):
	lowess_xy <- function(x, ..., var=c("Center.x", "Center.y"), xynames=c("Center.x_lowess", "Center.y_lowess")){
		# Smooth first to get the computed spar:
		fit_x <- lowess(x$UNIXtime, x[[var[1]]], ...)
		fit_y <- lowess(x$UNIXtime, x[[var[2]]], ...)
	
		y <- setNames(list(fit_x$y, fit_y$y), xynames)
		x <- cbind(x, y)
		return(x)
	}
	
	#Calculate the dot product and incidence angle
	getAngle <- function(x, y){
		if(length(dim(x))==0){
		  x <- t(x)
		}
		if(length(dim(y))==0){
		  y <- t(y)
		}
		lx <- sqrt(rowSums(x^2))
		ly <- sqrt(rowSums(y^2))
		acos( rowSums(x * y) / (lx * ly)) * 180/pi
	}
	
	dist2d <- function(x, p1, p2) {
		# Function for the distance between a point and a line between two points:
		dist2dOne <- function(x, p1, p2) {
			x <- unlist(x, use.names=FALSE)
			p1 <- unlist(p1, use.names=FALSE)
			p2 <- unlist(p2, use.names=FALSE)
			v1 <- p1 - p2
			v2 <- x - p1
			m <- cbind(v1,v2)
			abs(det(m))/sqrt(sum(v1*v1))
		} 
		#browser()
		#if(ncol(p1)==2){
		#	p2 <- p1[-1,]
		#	p1 <- p1[-nrow(p1),]
		#}
		if(nrow(p1) == nrow(x)-1){
			p1 <- rbind(p1, p1[nrow(p1),])
			p2 <- rbind(p2, p2[nrow(p2),])
		}
		if(nrow(p1) < nrow(x) || nrow(p1)==1){
			p1 <- matrix(p1, ncol=2, nrow=nrow(x), byrow=TRUE)
			p2 <- matrix(p2, ncol=2, nrow=nrow(x), byrow=TRUE)
		}
		sapply(seq_len(nrow(x)), function(i) dist2dOne(x[i,], p1[i,], p2[i,]))
	} 
	
	
	##### 0. Read the data: #####
	x <- read.table(x, sep="", header=TRUE, stringsAsFactors=FALSE, nrows=nrows)


	##### 1. Time: #####
	# Reshape time to R format. Number of seconds since 1970: #####
	x$UNIXtime <- unclass(as.POSIXct(paste(x$Date, x$Time), format="%Y-%m-%d %H:%M:%OS", tz="GMT"))
	# Convert to POSIX:
	x$DateTime<-as.POSIXct(x$UNIXtime, origin="1970-01-01", tz="GMT")
	# Sort the sonar data by ID and time
	x <- x[order(x$Id, x$UNIXtime), ]
	# Remove duplicates:
	if(rm.duplicates){
		x <- x[!duplicated(x$UNIXtime), ]
	}
	
	
	##### 2. X,Y: #####
	# Add Cartesian coordinates for the school and vessel:
	lat0 <- mean(x$Center.lat)
	lon0 <- mean(x$Center.lon) 
	schoolXY <- geoXY(x$Center.lat, x$Center.lon, lat0, lon0, unit=1) #distance (m) for each detection to mean postion
	x$Center.x <- schoolXY[,1]
	x$Center.y <- schoolXY[,2]
	shipXY <- geoXY(x$Ship.lat, x$Ship.lon, lat0, lon0, unit=1) #distance (m) for each detection to mean postion
	x$Ship.x <- shipXY[,1]
	x$Ship.y <- shipXY[,2]
	
	# Add the range from the vessel:
	x$Center.range <- sqrt( (x$Center.x - x$Ship.x)^2 + (x$Center.y - x$Ship.y)^2 )
	
	
	# Add the current information:
	if(length(currentSpeed)>0 && length(currentAngle)>0){
		# Convert to m/s:
		x$currentSpeed <- currentSpeed
		x$currentAngle <- currentAngle
	}
	else{
		x$currentSpeed <- NA
		x$currentAngle <- NA
	}
	
	
	##### 3. Filter: #####
	applyFilter <- function(x, filter){
		applyFilterOne <- function(var, x, filter){
			if(length(x[[var]]) && length(filter[[var]])==2){
				x[[var]] >= filter[[var]][1] & x[[var]] <= filter[[var]][2]
			}
			else{
				!logical(nrow(x))
			}
		}
		filter <- filter[names(filter) %in% names(x)]
		if(length(filter) && is.list(filter)){
			valid <- sapply(names(filter), applyFilterOne, x=x, filter=filter)
			valid <- apply(valid, 1, all)
			x <- x[valid, ]
		}
		
		return(x)
	}
	
	# Apply the filters before smoothing:
	x <- applyFilter(x, filter)
	#if(length(filter) && is.list(filter)){
	#	valid <- sapply(names(filter), applyFilter, filter=filter, x=x)
	#	valid <- apply(valid, 1, all)
	#	x <- x[valid, ]
	#}
	
	numPings <- table(x$Id)
	numPings <- rep(numPings, numPings)
	x$numPings <- numPings
	
	# Apply any numPings filter before smoothing:
	x <- applyFilter(x, filter)
	
	
	##### 3. Spline smooth: #####
	runSmoothing <- function(Id=NULL, data, sparfact, spar0){
		# Get the current school:
		if(length(Id)){
			data <- data[data$Id==Id, ]
		}
		
		# Smooth the school positions using spline:
		data <- smooth.spline_xy(data, sparfact=sparfact[1], spar0=spar0, xynames=c("Center.x_spline", "Center.y_spline"))
	
		# For the vessel, using only the default smoothing here:
		data <- smooth.spline_xy(data, sparfact=NULL, spar0=spar0, var=c("Ship.x", "Ship.y"), xynames=c("Ship.x_spline", "Ship.y_spline"))
	
	
		##### 4. Lowess smooth: #####
		data <- lowess_xy(data, xynames=c("Center.x_lowess", "Center.y_lowess"), ...)
	
		# For the vessel, using 10 values:
		f <- 10/nrow(data)
		data <- lowess_xy(data, var=c("Ship.x", "Ship.y"), xynames=c("Ship.x_lowess", "Ship.y_lowess"), f=f)
	
	
		##### 5. Heading and speed from diffing the smoothed positions: #####
		data <- getDiffHeadingSpeed(data, var="Center", type="")
		data <- getDiffHeadingSpeed(data, var="Ship", type="")
	
		data <- getDiffHeadingSpeed(data, var="Center", type="spline")
		data <- getDiffHeadingSpeed(data, var="Ship", type="spline")
	
		data <- getDiffHeadingSpeed(data, var="Center", type="lowess")
		
		
		seq1 <- c(seq_len(nrow(data)-1), nrow(data)-1)
		seq2 <- seq1 + 1
		d <- dist2d(data[c("Center.x", "Center.y")], p1=data[seq1, c("Ship.x_spline", "Ship.y_spline")], p2=data[seq2, c("Ship.x_spline", "Ship.y_spline")])
		#data$DistToTrack <- median(d)
		data$DistToTrack <- d
		
		
		
		return(data)
	}
	
	if(smooth){
		if(byId){
			uId <- unique(x$Id)
			temp <- TSD::papply(uId, runSmoothing, data=x, sparfact=sparfact, spar0=spar0, cores=cores)
			#x <- do.call(rbind, temp)
			x <- as.data.frame(data.table::rbindlist(temp))
		}
		else{
			x <- runSmoothing(data=x, sparfact=sparfact, spar0=spar0)
		}
	}
	else{
		return(x)
	}
	
	# Apply the filters after smoothing:
	x <- applyFilter(x, filter)
	
	
	##### 6. Incidence angle: #####
	# Get the vector from the school to the vessel:
	x$diffAng <- getAngle(cbind(x$Center.dx_lowess, x$Center.dy_lowess), cbind(x$Ship.dx_spline, x$Ship.dy_spline))
	x$incAng_lowess <- getAngle(cbind(x$Center.dx_lowess, x$Center.dy_lowess), cbind(x$Ship.x_spline - x$Center.x_lowess, x$Ship.y_spline - x$Center.y_lowess))
	x$incAng_spline <- getAngle(cbind(x$Center.dx_spline, x$Center.dy_spline), cbind(x$Ship.x_spline - x$Center.x_spline, x$Ship.y_spline - x$Center.y_spline))
	
	
	##### 7. Rotation angle: #####
	r <- rotate3D(cbind(x=x$Ship.x_spline - x$Center.x_lowess, y=x$Ship.y_spline - x$Center.y_lowess, z=0), "z", x$Center.heading_lowess, paired=TRUE)
	r <- atan2(r[,2], r[,1])
	
	x$ShipAng_lowess <- r
	
	x$incRot_lowess <- ang2rot(r)
	
	if(tail(x$incRot_lowess, 1) < head(x$incRot_lowess, 1)){
		warning("Incidence angle rotation multiplied by -1 to obtain time along the positive x axis")
		x$incRot_lowess <- -1 * x$incRot_lowess
	}
	
	return(x)
}