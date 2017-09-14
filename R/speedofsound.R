#*********************************************
#*********************************************
#' Computes the speed of sound of the ocean for the given temperature, salinity and pressure values. Method adopted from "sw_svel.m" and based on UNESCO 1983.
#'
#' @param T  is the temperature in degrees celcius, or a list of the inputs in the TSD file format.
#' @param S  are the salinity values in parts per thousand.
#' @param P  are the pressure values in decibar or Pascal.
#' @param Z  are the z-positions in meters (negative). Must be given along with 'P0', 'rho' and 'g'.
#' @param P0  is the air pressure in decibar or Pascal (see 'Pain').
#' @param rho  is the mass density of the sea in kg pr m^3.
#' @param g  is the gravitational acceleration in meters per s^2.
#' @param Pain  is TRUE if pressure is given in Pascal and FALSE if given in decibar (10000 Pascal).
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname speedofsound
#'
speedofsound<-function(T, S, P, Z=NULL, P0=10.1325, rho=1027, g=9.82, Pain=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2008-02-23 - Finished.
	# Update: 2009-06-12 - Changed to support list and matrix input.
	# Update: 2010-02-02 - Corrected bug when specifying pressure, and added options for specifying the air pressure and mass density of the sea water. Added TSD input support.
	# Last: 2010-02-07 - Method adopted from "sw_svel.m". Old version renamed to speedofsound_crude and moved to the "-unused" directory.
	########### DESCRIPTION: ###########
	# Computes the speed of sound of the ocean for the given temperature, salinity and pressure values. Method adopted from "sw_svel.m" and based on UNESCO 1983.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---T--- is the temperature in degrees celcius, or a list of the inputs in the TSD file format.
	# ---S--- are the salinity values in parts per thousand.
	# ---P--- are the pressure values in decibar or Pascal.
	# ---Z--- are the z-positions in meters (negative). Must be given along with 'P0', 'rho' and 'g'.
	# ---P0--- is the air pressure in decibar or Pascal (see 'Pain').
	# ---rho--- is the mass density of the sea in kg pr m^3.
	# ---g--- is the gravitational acceleration in meters per s^2.
	# ---Pain--- is TRUE if pressure is given in Pascal and FALSE if given in decibar (10000 Pascal).
	
	
	##################################################
	##################################################
	##### Preparation #####
	if(is.list(T)){
		# If speed is present in the input list as the TSD-variable 'isps', this is simply returned:
		if(!is.null(T$isps)){
			return(T$isps)
			}
		# Else names of the input list are used to identify the reqiured variables:
		namesT=names(T)=tolower(names(T))
		# Salinity:
		namesS=match(c("slty","s"),namesT)
		S=T[[namesS[1]]]
		# Pressure:
		namesP=match(c("ihpr","p"),namesT)
		P=T[[namesP[1]]]
		# Depth:
		namesZ=match(c("pszc","z"),namesT)
		Z=T[[namesZ[1]]]
		# Mass density:
		namesrho=match(c("rho0","rho"),namesT)
		rho=T[[namesrho[1]]]
		# Gravitational acceleration
		namesg=match(c("gacc","g"),namesT)
		g=T[[namesg[1]]]
		# Pressure at sea level:
		namesP0=match(c("hpr0","p0","airpr"),namesT)
		P0=T[[namesP0[1]]]
		# Temperature:
		namesT=match(c("temp","t"),namesT)
		T=T[[namesT[1]]]
		}
	# Calculating the pressure 'P' if needed:
	if(is.null(P) && !any(is.null(Z),is.null(rho),is.null(g),is.null(P0))){
		if(Pain){
			P=P0-rho*g*Z
			}
		else{
			P=P0-rho*g*Z*1e-4
			}
		}
	# convert 'P' to bars as used in UNESCO routines
	if(Pain){
		P = P*1e-5
		}
	else{
		P = P/10
		}
	T68 = T * 1.00024
	
	
	##### Execution #####
	#------------
	# eqn 34 p.46
	#------------
	c00 =  1402.388
	c01 =  5.03711
	c02 = -5.80852e-2
	c03 =  3.3420e-4
	c04 = -1.47800e-6
	c05 =  3.1464e-9
	
	c10 =  0.153563
	c11 =  6.8982e-4
	c12 = -8.1788e-6
	c13 =  1.3621e-7
	c14 = -6.1185e-10
	
	c20 =  3.1260e-5
	c21 = -1.7107e-6
	c22 =  2.5974e-8
	c23 = -2.5335e-10
	c24 =  1.0405e-12
	
	c30 = -9.7729e-9
	c31 =  3.8504e-10
	c32 = -2.3643e-12
	
	Cw = ((((c32*T68 + c31)*T68 + c30)*P + ((((c24*T68 + c23)*T68 + c22)*T68 + c21)*T68 + c20))*P + ((((c14*T68 + c13)*T68 + c12)*T68 + c11)*T68 + c10))*P + ((((c05*T68 + c04)*T68 + c03)*T68 + c02)*T68 + c01)*T68 + c00
	
	#-------------
	# eqn 35. p.47
	#-------------
	a00 =  1.389
	a01 = -1.262e-2
	a02 =  7.164e-5
	a03 =  2.006e-6
	a04 = -3.21e-8
	
	a10 =  9.4742e-5
	a11 = -1.2580e-5
	a12 = -6.4885e-8
	a13 =  1.0507e-8
	a14 = -2.0122e-10
	
	a20 = -3.9064e-7
	a21 =  9.1041e-9
	a22 = -1.6002e-10
	a23 =  7.988e-12
	
	a30 =  1.100e-10
	a31 =  6.649e-12
	a32 = -3.389e-13
	
	A = ((((a32*T68 + a31)*T68 + a30)*P + (((a23*T68 + a22)*T68 + a21)*T68 + a20))*P + ((((a14*T68 + a13)*T68 + a12)*T68 + a11)*T68 + a10))*P + (((a04*T68 + a03)*T68 + a02)*T68 + a01)*T68 + a00
	
	#------------
	# eqn 36 p.47
	#------------
	b00 = -1.922e-2
	b01 = -4.42e-5
	b10 =  7.3637e-5
	b11 =  1.7945e-7
	
	B = b00 + b01*T68 + (b10 + b11*T68)*P
	
	#------------
	# eqn 37 p.47
	#------------
	d00 =  1.727e-3
	d10 = -7.9836e-6
	
	D = d00 + d10*P
	
		
	##### Output #####
	#------------
	# eqn 33 p.46
	#------------
	Cw + A*S + B*S*sqrt(S) + D*S^2
	##################################################
	##################################################
	}
