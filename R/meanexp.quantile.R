#*********************************************
#*********************************************
#' Estimates the parameter of an exponentially distributed 'X' based on the quantile of values. F(x) = 1 - exp(-x/beta) = p => betahat = quantile(x,p) / log(1/(1-p)). Arrays up to 3 dimensions are accepted.
#'
#' @param x  is the data.
#' @param prob  is the probability.
#' @param type  is the type of the quantile calculation (see quantile()).
#' @param MARGIN  is used in the same way as in apply().
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom stats quantile
#'
#' @export
#' @rdname meanexp.quantile
#'
meanexp.quantile=function(x,prob=0.5,type=6,MARGIN=NULL){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-05-23 - Clean version.
	########### DESCRIPTION: ###########
	# Estimates the parameter of an exponentially distributed 'X' based on the quantile of values. F(x) = 1 - exp(-x/beta) = p => betahat = quantile(x,p) / log(1/(1-p)). Arrays up to 3 dimensions are accepted.
	# The type=6 is chosen to minimize the bias, and seems suited for the exponential distribution. This function replaces meanexp.conditional() in the noise estimation in echoIBM.
	# Takes one hour:
	# f=function(n=1e4,seed=NULL){
	# 	m=zeros(9,n)
	# 	if(length(seed)==0){
	# 		seed=round(runif(1)*1e6)
	# 		}
	# 	 set.seed(seed)
	# 	 for(r in seq_len(n)){
	# 	 	if((r%%1000) == 0){
	# 			cat(r,"\t")
	# 	 		}
	# 	 	for(i in seq_len(9)){
	# 	 		x=rexp(1e3)
	# 	 		m[i,r]=meanexp.quantile(x,prob=0.2,type=i)
	# 	 		}
	# 	 	}
	# 	 abs(rowMeans(m)-1)
	# 	}
	# n=10
	# m=zeros(9,n)
	# seeds=round(runif(10)*1e6)
	# seeds
	# 
	# for(i in seq_len(n)){
	# 	cat("Run ",i,"\n")
	# 	system.time(m[,i]<-f(1e5,seed=seeds[i]))
	# 	}
	# ploto(m[,1])
	# for(i in 2:n){
	# 	lines(m[,i],col=i,type="o")
	#	}
	# lines(rowMeans(m),lwd=3,type="o")
	# lines(apply(m,1,median),lwd=3,type="o",col="green")
	# This shows that type=1, 3, 4, and 6 are good, so we choose 6, which is sophisitcated and used by SPSS.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is the data.
	# ---prob--- is the probability.
	# ---type--- is the type of the quantile calculation (see quantile()).
	# ---MARGIN--- is used in the same way as in apply().
	
	
	##################################################
	##################################################
	########## Preparation ##########
	# Treat vectors:
	if(length(dim(x))==0){
		return(quantile(x,prob,type=type,na.rm=TRUE) / log(1/(1-prob)))
		}
	# Treat matrices:
	if(length(dim(x))==2 && length(MARGIN)==0){
		MARGIN=2
		}
	if(length(dim(x))==3 && length(MARGIN)==0){
		MARGIN=2:3
		}
	
	
	########## Execution and output ##########
	# Output:
	apply(x,MARGIN,function(x) quantile(x,prob,type=type,na.rm=TRUE)) / log(1/(1-prob))
	##################################################
	##################################################
	}
