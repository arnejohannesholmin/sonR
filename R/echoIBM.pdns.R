#*********************************************
#*********************************************
#' Calculates the periodic noise for all voxels based on estimates of the parameters of the periodic noise for ONE SINGLE TIME STEP. The function assumes that the phase is correctly given in data$pns3.
#'
#' @param data  is a list containing beam configuration, and periodic noise parameters. Speficically the following variables must be included: Beams variables: 'sint', 'lenb'; periodic noise variables: 'pns1', 'pns2', 'pns3', 'acfq', 'harm'; and either 'bgns' or 'numb'.
#' @param indt  is the time step identifier, which if given selects the phase parameter from the matrix data$pn3M.
#' @param TVG  is TRUE if TVG should be added to the periodic noise.
#' @param TVG.exp  is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw apply.TVG
#' @importFrom TSD zeros
#'
#' @export
#' @rdname echoIBM.pdns
#'
echoIBM.pdns<-function(data, indt=NULL, TVG=FALSE, TVG.exp=2){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2012-06-26 - Clean version.
	########### DESCRIPTION: ###########
	# Calculates the periodic noise for all voxels based on estimates of the parameters of the periodic noise for ONE SINGLE TIME STEP. The function assumes that the phase is correctly given in data$pns3.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---data--- is a list containing beam configuration, and periodic noise parameters. Speficically the following variables must be included: Beams variables: 'sint', 'lenb'; periodic noise variables: 'pns1', 'pns2', 'pns3', 'acfq', 'harm'; and either 'bgns' or 'numb'.
	# ---indt--- is the time step identifier, which if given selects the phase parameter from the matrix data$pn3M.
	# ---TVG--- is TRUE if TVG should be added to the periodic noise.
	# ---TVG.exp--- is the exponent of the eamotric spreading of the sound wave, theoretically 2 for Sv and 4 for TS.
	
	
	##################################################
	##################################################
	########## Preparation ##########
	if(length(dim(data$bgns))==2){
		data$pdns=zeros(dim(data$bgns))
		}
	else if(length(data$lenb)>0 && length(data$numb)>0){
		data$pdns=zeros(max(data$lenb),data$numb)
		}
	else{
		data$pdns=0
		}
	data$pnsA=data$pdns
	
	
	########## Execution ##########
	if(!all(c("sint","lenb","numb","bgns","acfq","badb","harm","pns1","pns2") %in% names(data)) && !any(c("pn3M","pns3") %in% names(data))){
		warning("Periodic noise not subtracted. Specificaiton of periodic noise not present in the data")
		return(data$pdns)
		}
	periodic=which(data$badb==1)
	# Define the frequency of the periodic noise, and the sequence along beams:
	a=2*pi * data$acfq*data$sint * data$harm[periodic]
	r=seq_len(max(data$lenb))
	# Calculate the periodic noise:
	if(length(indt)>0){
		if(indt[1]<1 && indt[1]>NCOL(data$pn3M)){
			warning(paste("'indt' (",indt[1],") chosen outside of the valid range [1,",NCOL(data$pn3M),"]",sep=""))
			}
		data$pns3=data$pn3M[,indt[1]]
		}
	data$pdns[,periodic] = t( data$pns1[periodic] * 10^( data$pns2[periodic] * sin( outer(a,r) + data$pns3[periodic] ) ) )
	data$pnsA[,periodic] = matrix(data$pns1[periodic]*10^(data$pns2[periodic]),nrow=max(data$lenb),ncol=length(periodic),byrow=TRUE)
			
	# Add TVG if required:
	if(TVG){
		data$pdns[,periodic] = apply.TVG(data$pdns, data, TVG.exp=TVG.exp)[,periodic]
		data$pnsA[,periodic] = apply.TVG(data$pnsA, data, TVG.exp=TVG.exp)[,periodic]
		}
	
	
	########## Output ##########
	data[c("pdns","pnsA")]
	##################################################
	##################################################
	}
