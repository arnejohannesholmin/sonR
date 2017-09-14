#*********************************************
#*********************************************
#' Parses segmentaiton file names and extracts segmentation parameter names and values.
#'
#' @param x  is a file name containing segmentation parameter names and values in pairs enclosed in brackets <> (example: "S2009116_D20091113_E0001_<bwGp><1><lsth><1e-05><usth><2e-04><sgth><1.0e-01>_T081147.seg").
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD even NAs odd
#'
#' @export
#' @rdname segfilename2segpar
#'
segfilename2segpar<-function(x,force.numeric=TRUE){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2013-07-16 - Clean version.
	########### DESCRIPTION: ###########
	# Parses segmentaiton file names and extracts segmentation parameter names and values.
	########## DEPENDENCIES: ###########
	#
	############ VARIABLES: ############
	# ---x--- is a file name containing segmentation parameter names and values in pairs enclosed in brackets <> (example: "S2009116_D20091113_E0001_<bwGp><1><lsth><1e-05><usth><2e-04><sgth><1.0e-01>_T081147.seg").
	
	
	##################################################
	##################################################
	##### Preparation #####
	# Locate the crocodiles:
	cl=gregexpr("\\{",x)[[1]]
	cr=gregexpr("\\}",x)[[1]]
	
	# Get the positions of valied separators, so that everything between a "<" and a ">" is read as one input:
	this=0
	croc=cl[1]
	iter=1
	while(this<max(cl,cr)){
		if(even(iter)){
			this=min(cl[cl>croc[iter]])
			}
		else{
			this=min(cr[cr>croc[iter]])
			}
		iter=iter+1
		croc=c(croc,this)
		}
	croc=matrix(croc,byrow=TRUE,ncol=2)
		
	# If the number of pairs is still odd, issue a warning:
	N=nrow(croc)
	if(odd(N)){
		warning("The file name should contain paris of variable names and values such as \"<code>\"\"<5>\" ")
		}
	
	
	##### Execution and output #####
	data=NAs(N)
	Nseq=seq_len(N)
	# Extract variable names and values:
	for(i in Nseq){
		data[i]=substr(x,croc[i,1]+1,croc[i,2]-1)
		}
	
	# Create the output list:
	datanames=data[odd(Nseq)]
	# Find the numeric values:
	arenumeric = !is.na(as.numeric(data[even(Nseq)]))
	# Convert to a list and convert to numeric:
	if(force.numeric){
		data=as.list(as.numeric(data[even(Nseq)]))
		}
	else{
		data=as.list(data[even(Nseq)])
		data[arenumeric]=lapply(data[arenumeric],as.numeric)
		}
	names(data)=datanames[seq_along(data)]
	
	# Split 'bwGp' into an x, y, and z component:
	bwGp = names(data)=="bwGp"
	if(length(bwGp)>0){
		data$bwGp_x=data$bwGp
		data$bwGp_y=data$bwGp
		data$bwGp_z=data$bwGp
		}
	else if(any(length(data$bwGp_x)>0, length(data$bwGp_y)>0, length(data$bwGp_z)>0)){
		data$bwGp=mean(unlist(data[c("bwGp_x","bwGp_y","bwGp_z")]))
		}
	
	
	# Return:
	data
	##################################################
	##################################################
	}
