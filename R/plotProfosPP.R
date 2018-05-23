#*********************************************
#*********************************************
#' Plot output from readProfosPP.
#'
#' This function creates two plots, (1) one showing the raw and spline smoothed track of the vessel (ship) and the raw and spline and lowess smoothed track of the school, and (2) a related plot showing the mean volume backscattering coefficient/strength of the school as a function of incidence angle rotation (increasing along the unit circle). 
#'
#' @param x					The output from \code{readProfosPP}.
#' @param ind				A numeric vector indicating which plots to produce, 1 being the track plot and 2 the sv plot.
#' @param col,pch,cex		A list of color, ploting character and character expansion values for the ship, school and links between these at transitions between weak and strong echo.
#' @param lwd				A list of line width values for the (smoothed) tracks (plotted with type="o"), the links, and the sv.
#' @param plotlink			Logical: If TRUE plot links.
#' @param linkvar			The type of smoothing of the school to connect ship positions to, one of "spline" (less smoothing) and "lowess" (stronger smoothing).
#' @param every				Plot links at every \code{every} ping.
#' @param use				Whether to use linear ('sv') or logarithmical ('Sv') mean volume backscattering.
#' @param svQ				A vector of length 2 defining the quantiles used for defining the amplitude of the sine wave which is fitted to the mean Sv to obtain transitions between weak and strong echo.
#' @param ...				Arguments passed to \code{\link{plot}}.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname readProfosPP
#'
plotProfosPP <- function(x, 
	ind = c(1, 2), 
	col = list(ship=c("black", "purple"), school=c("black", "brown", "orange"), link="green4"), 
	pch = list(ship="*", school=1), 
	cex = list(ship=c(1, 0.5), school=c(1, 0.5, 0.5), link=1), 
	lwd = list(track=0.5, link=1, sv=1), 
	plotlink = TRUE, 
	linkvar = "lowess", 
	every = NULL, 
	use = c("sv", "Sv"), 
	svQ = c(0.1, 0.9),
	...){
	
	# Define color blind safe palette (not used any more):
	#cbc <- c(0, 1), col=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	lll <- list(...)

	
	repSpec <- function(x){
		x$ship <- rep(x$ship, length.out=2)
		x$school <- rep(x$school, length.out=3)
		x
	}
	
	col <- repSpec(col)
	pch <- repSpec(pch)
	cex <- repSpec(cex)
	
	sv <- x$Sv.mean
	if(use[1]=="sv"){
		sv <- 10^(sv/10)
	}
		
	time <- as.POSIXct(x$UNIXtime, origin="1970-01-01", tx="UTC")
	
	# Fit a sine wave to the logarithmic Sv:
	sin_phase <- function(phase, x){
		sin( 2 * (x$incRot_lowess - phase))
	}
	# The funciton to optimize:
	f <- function(par, x, lev, fact){
		phase <- par[1]
		x1 <- (x$Sv.mean - lev) / fact
		#x2 <- sin_a_phase(a=a, phase=phase, x=x)
		x2 <- sin_phase(phase=phase, x=x)
		sum((x1-x2)^2)
	}

	# Fit a sine wave to the mean Sv as a function of rotation angle:
	lev <- mean(x$Sv.mean)
	fact <- diff(quantile(x$Sv.mean, svQ)) / 2
	o <- optimize(f, c(-pi, pi), x=x, lev=lev, fact=fact)
	# Get the fitted phase and create a vector at each intersection with 0 in the size wave, corresponding to transitions between weak and strong echo in the shcool (45, 135, 225, 315 degrees)
	start <- o$minimum
	s <- seq(floor(min(x$incRot_lowess) / (pi/2)), ceiling(max(x$incRot_lowess) / (pi/2)))
	v <- start + s * pi / 2
	inside <- which(v > min(x$incRot_lowess) & v < max(x$incRot_lowess))
	
	# Find the corresponding pings:
	inds <- apply(abs(outer(x$incRot_lowess, v, "-")), 2, which.min)
	# Discard values outside of the range of rotation angles for the labels:
	inds <- inds[inside]
	linkText <- seq_along(inside)
	

	# Plot the smoothing of the school
	if(1 %in% ind){
		if(is.list(lll$xlim)){
			thisxlim <- xlim[[1]]
		}
		else{
			thisxlim <- range(x$Center.x, x$Center.x_spline, x$Ship.x)
		}
		if(is.list(lll$ylim)){
			thisylim <- ylim[[1]]
		}
		else{
			thisylim <- range(x$Center.y, x$Center.y_spline, x$Ship.y)
		}
		
		
		#Id123 <- match( x$Id, unique(x$Id))
		
		# Vessel:
		plot(x$Ship.x, x$Ship.y, xlim=thisxlim, ylim=thisylim, xlab="x", ylab="y", col=col$ship[1], pch=pch$ship[1], cex=cex$ship[1], ...)
		points(x$Ship.x_spline, x$Ship.y_spline, type="o", lwd=lwd$track, col=col$ship[2], pch=pch$ship[2], cex=cex$ship[2])
		
		if(length(every)){
			sx1 <- x[[paste0("Center.x_", linkvar)]]
			sy1 <- x[[paste0("Center.y_", linkvar)]]
			valid <- c(seq(1, length(sx1), every), length(sx1))
			segments(sx1[valid], sy1[valid], x$Ship.x_spline[valid], x$Ship.y_spline[valid], lwd=lwd$link)
		}
		
		# School:
		points(x$Center.x, x$Center.y, col=col$school[1], pch=pch$school[1], cex=cex$school[1])
		lines(x$Center.x_spline, x$Center.y_spline, type="o", lwd=lwd$track, col=col$school[2], pch=pch$school[2], cex=cex$school[2])
		lines(x$Center.x_lowess, x$Center.y_lowess, type="o", lwd=lwd$track, col=col$school[3], pch=pch$school[3], cex=cex$school[3])
		
		# Add lines between the school and the vessel corresponding to the vertical lines in the mean sv plot:
		if(plotlink){
			sx1 <- x[[paste0("Center.x_", linkvar)]][inds]
			sy1 <- x[[paste0("Center.y_", linkvar)]][inds]
			sx2 <- x$Ship.x_spline[inds]
			sy2 <- x$Ship.y_spline[inds]
			segments(sx1, sy1, sx2, sy2, lwd=lwd$link, col=col$link)
			points(c(sx1, sx2), c(sy1, sy2), pch=20, col=col$link)
			towardsShip <- 0.8
			tx <- sx1 + (sx2 - sx1) * towardsShip
			ty <- sy1 + (sy2 - sy1) * towardsShip
			text(tx, ty, linkText, col=col$link, cex=cex$link)
		}
	}
	
	
	#if(2 %in% ind){
	#	plot(time, x$Center.heading * 180/pi, type="o", col=1, lwd=0.5)
	#	lines(time, x$Center.heading_lowess * 180/pi, type="o", col=1, lwd=0.5)
	#	lines(time, x$Center.heading_spline * 180/pi, type="o", col=3, lwd=0.5)
	#}
	#
	#if(3 %in% ind){
	#	plot(time, x$Sv.mean, type="o", col=1, lwd=0.5)
	#}
	
	#if(2 %in% ind){
	#	plot(x$incAng_spline, type="o")
	#	lines(x$incAng_lowess, type="o", col=3, lwd=0.5)
	#}
	#if(3 %in% ind){
	#	plot(time, x$incRot_lowess, type="o", col=1, lwd=0.5)
	#}
	
	
	if(2 %in% ind){
		ang <- x$incRot_lowess * 180/pi
		
		if(is.list(lll$xlim)){
			thisxlim <- xlim[[2]]
		}
		else{
			thisxlim <- range(ang)
		}
		if(is.list(lll$ylim)){
			thisylim <- ylim[[2]]
		}
		else{
			thisylim <- range(sv)
		}
		
		plot(ang, sv, col=1, lwd=lwd$sv, xlim=thisxlim, ylim=thisylim, xlab="Rotation angle", ylab=paste("Mean", use[1]), type="l", xaxt="n", ...)
		axis(1, at=s * 90, s * 90, las=2)
		
		# Add vertical lines at the maxima:
		v <- v * 180/pi
		abline(v = v, col=col$link)
		if(plotlink){
			text(v[inside], thisylim[2] - diff(thisylim) * 0.05, linkText, col=col$link, cex=cex$link)
		}
	}
	
	return(NULL)
}
