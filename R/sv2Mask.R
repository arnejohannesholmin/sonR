#' Bayesian segmentation method for one ping
#' 
#' Function to get the segmentation mask of one ping from the Bayesian segmentation method (Holmin et al., 2024)
#' 
#' @param sv  A matrix of sv (linear values) for one ping.
#' @param betaN  A vector of background noise for all beams.
#' @param beta0  The lower schooling threshold given in linear value. Suggested value: 1e-5.
#' @param beta1  The upper schooling threshold given in linear value. Suggested value: 1e-2.
#' @param h  The standard deviation in the Gaussian kernel used to smooth the probability of "not-school". Suggested value: 5.
#' @param c  The segmentation threshold. Suggested value: 1e-5.
#' @param beams  A list of beam configuration data. Must contain he following elements: lenb, numb, plsl, asps, sint, absr, rofs. Read this from aa .beams file in the TSD format.
#'
#' @importFrom SimradRaw apply.TVG
#' @importFrom expint expint_E1
#'
#' @export
#'
sv2MaskOnePing <- function(sv, betaN, beta0, beta1, h, c, beams) {

	# Expand the lower schooling threshold beta0, the upper schooling threshold beta1, segmentation threshold c, and the noise betaN to a matrix with dimensions (number of samples along beams, number of beams):
	beta0 <- matrix(beta0, nrow = beams$lenb, ncol = beams$numb, byrow = TRUE)
	beta1 <- matrix(beta1, nrow = beams$lenb, ncol = beams$numb, byrow = TRUE)
	c <- matrix(c, nrow = beams$lenb, ncol = beams$numb, byrow = TRUE)
	
	betaN <- matrix(betaN, nrow = beams$lenb, ncol = beams$numb, byrow = TRUE)
	
	# Add TVG using the SimradRaw package:
	betaNWithTVG <- SimradRaw::apply.TVG(betaN, beams, TVG.exp = 2)
	
	# Impose the restriction that beta0 must be larger or equal to the background noise:
	beta0BelowNoise <- beta0 < betaNWithTVG
	beta0[beta0BelowNoise] <- betaNWithTVG[beta0BelowNoise]	
	


	# Run the cummulative distribution function evaluated at beta0, resulting in the probability of "not school":
	p <- probabilityOfNotSchool(
		beta0 = beta0,
		beta1 = beta1,
		sv = sv, 
		betaN = betaN
	)
	
	# Smooth in the logarithmic domain:
	pSmooth <- exp(apply(log(p), 2, GaussianSmoothOneDimension, h = h))
	
	# Get the segmentation mask by thresholding the smoothed probabilities of "not school":
	mask <- pSmooth < c

	
	# Return the output:
	return(mask)
}

# Define the function for calculating the cummulative distribution. This function takes up to one second for one ping of MS70 data, mostly due to the expint_E1() function:
probabilityOfNotSchool <- function(beta0, beta1, sv, betaN){
	
	# The ratio of exponential integrals is Eq. (6) in "Bayesian segmentation method.pdf":
	H <- function(b, sv, betaN, beta1){
		expint::expint_E1(sv / (b + betaN)) / expint::expint_E1(sv / (beta1 + betaN))
	}
	
	# Make the ratio:
	out <- H(beta0, sv = sv, betaN = betaN, beta1 = beta1) / 
		H(beta1, sv = sv, betaN = betaN, beta1 = beta1)
	
	return(out)
}

# Function for Gaussian smoothing along a vector:
GaussianSmoothOneDimension <- function(x, h, nh = 3) {
	lag <- seq(nh * -h, nh * h)
	GaussianKernel <- dnorm(lag, 0, h)
	filter(x, filter = GaussianKernel)
}


