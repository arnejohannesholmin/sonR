#*********************************************
#*********************************************
#' Estimates the parameter 'beta' of an exponentially distributed 'X' based on the quantile of values. 
#'
#' This function accepts an addition to the expected value of 'X', 'a', which must be known. The calculation is an optimizatoin of 'beta' from an expression deduced from the cummulative distribution function (CDF) of an arbitrary 'X' (which is the mean of the CDF of X_t):
#' X_t ~ exp(beta_t), where beta_t = beta + a_t
#' F_X(x) = 1/N * sum_t=1^N(F_X_t(x_t)) = 1 - 1/N * sum_t=1^N(exp(-x/beta_t)),
#' which is estimated by
#' u = 1 - 1/N * sum_t=1^N(exp(-Q(u)/beta_t))
#' u = 1 - 1/N * sum_t=1^N(exp(-Q(u)/(beta+a_t))), where Q(u) is the u-quantile.
#' This expression is solved numerically for 'beta'
#'
#' @param x				The data.
#' @param lower,upper	Used in \code{\link{optimize}}.
#' @param prob			The probability.
#' @param a				The added constant.
#' @param type			The type of the quantile calculation (see quantile()).
#' @param MARGIN		Used in the same way as in apply().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD ones zeros
#' @importFrom stats optimize quantile
#'
#' @export
#' @rdname meanexp.quantile_a
#'
meanexp.quantile_a <- function(x, lower=0, upper=1, prob=0.5, a=1, type=6, MARGIN=NULL){
	
	############### LOG: ###############
	# Start: 2012-05-23 - Clean version.
	# Last: 2012-09-23 - Changed to return a list of the estimates and the function values.
	
	# The expression which is the Exponential cummulative distribution function estimated by the quantile q(prob)
	CDF <- function(beta, q, a, prob, N){
		#abs(sum(exp(-q/(beta+a)))-N*(1-prob))
		# Change implemented after comment from Dag Tjostheim 2012-11-02:
		(sum(exp(-q/(beta+a)))-N*(1-prob))^2
	}
	getbeta <- function(x, lower, upper, prob, type, a){
		N <- length(x)
		Q <- quantile(x, prob, type=type, na.rm=TRUE)
		o <- optimize(CDF, c(lower,upper), q=Q, a=a, prob=prob, N=N)
		c(o$minimum, o$objective)
	}

	# Treat vectors:
	if(length(dim(x))==0){
		dim(x) <- c(length(x),1)
		MARGIN <- 2
	}
	# Treat matrices:
	if(length(dim(x))==2 && length(MARGIN)==0){
		MARGIN <- 2
	}
	if(length(dim(x))==3 && length(MARGIN)==0){
		MARGIN <- 2:3
	}
	
	# Output:
	out <- apply(x, MARGIN, getbeta, lower=lower, upper=upper, prob=prob, type=type, a=a)
	out <- split(out, c(cbind(zeros(length(out)/2),ones(length(out)/2))))
	list(minimum=out[[1]], objective=out[[2]])
}
