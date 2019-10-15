#*********************************************
#*********************************************
#' Computes the z-position based on CTD data (if missing).
#'
#' @param ctd  is the conductivity, temperature and depth data.
#' @param Pain  is TRUE if pressure "ihpr" is given in Pascal and FALSE if given in decibar relative to surface pressure (10000 Pascal, giving values approximately equivalent to water depth).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname getzfromctd
#'
getzfromctd<-function(ctd, Pain=FALSE){
	
	namesctd=names(ctd)=tolower(names(ctd))
	
	namesP=match(c("ihpr","p"),namesctd)
	P=ctd[[namesP[1]]]
	
	namesZ=match(c("pszc","z"),namesctd)
	Z=ctd[[namesZ[1]]]
	
	namesrho=match(c("rho0","rho"),namesctd)
	rho=ctd[[namesrho[1]]]
	
	namesg=match(c("gacc","g"),namesctd)
	g=ctd[[namesg[1]]]
	
	namesP0=match(c("hpr0","p0","airpr"),namesctd)
	P0=ctd[[namesP0[1]]]
	
	
	if(is.null(Z) && !any(is.null(P),is.null(rho),is.null(g),is.null(P0))){
		if(Pain){
			Z=(P0-P)/(rho*g)
		}
		else{
			Z=(P0-P)/(rho*g)*1e4
		}
	}
	Z
}
