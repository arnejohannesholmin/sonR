#*********************************************
#*********************************************
#' Extracts a subset of TSD data according to the ranges given for the variables of 'data' in 'range'. 
#'
#' @param data  is a list of elements named according to the TSD file format.
#' @param range  is a list of elements with names matching names if 'data', specifying the range of the corresponding elements.
#' @param this  is an optional previously generated list of indexes as returned from extract().
#' @param treated  A vector of variable names already treated.
#' @param ineq  is a string giving the inequality function to apply to the selection ("<" or "<=").
#' @param ind.out  is TRUE if the indexes for the elements segmented are to be returned.
#' @param insert.NA  is TRUE if the discarded data are to be kept as NA.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname extract.range.TSD
#'
extract.range.TSD <- function(data=list(), range=list(), this=list(), treated=NULL, ineq="<", ind.out=FALSE, insert.NA=FALSE){

	############### LOG: ###############
	# Start: 2013-06-06 - Clean version.
	
	##### Preparation #####
	# Get the names of 'range' and separate out the names concerinig voxels ('namesposx') and fish ('namesposf'):
	# If the mean volume backscattering strength (Sv) is to be subsetted, transform to volume backscattering coefficient (linear, sv):
	if("mvbs" %in% names(range)){
		range$vbsc=10^(range$mvbs/10)
		range$mvbs=NULL
		}
	namesrange=intersect(names(range),names(data))
	namesposx=intersect(namesrange,c("psxx","psyx","pszx"))
	namesposf=intersect(namesrange,c("psxf","psyf","pszf"))
	namesacoustic=intersect(namesrange,c("vbsc"))
	namesrange=setdiff(namesrange,c(namesposx,namesposf,namesacoustic))
	
		
	##### Execution #####
	# 2018-11-17: This seems like repeated code which sould be treated by a function:
	
	# Positions of the voxels, simultaneously subsetted:
	if(length(namesposx)>0){
		thisind=TRUE
		for(i in seq_along(namesposx)){
			thisind=thisind & do.call(ineq,list(min(range[[i]]),data[[namesposx[i]]])) & do.call(ineq,list(data[[namesposx[i]]],max(range[[i]])))
			}
		presentnames=intersect(c("psxx","psyx","pszx","vbsc","mvbs","volx","pr0s"),names(data))
		# Insert NA or subset the data:
		if(insert.NA){
			data[presentnames]=lapply(data[presentnames],function(x) {x[!thisind]=NA; x})
			}
		else{
			data[presentnames]=lapply(data[presentnames],function(x) x[thisind])
			}
		# Add to 'this$ind' (which are indices, not logical):
		this$ind[intersect(presentnames,treated)]=lapply(this$ind[intersect(presentnames,treated)],function(x) x[thisind])
		this$ind[setdiff(presentnames,treated)]=rep(list(which(thisind)),length(setdiff(presentnames,treated)))
		treated=union(treated,presentnames)
		}	
	
	# Positions of the fish, simultaneously subsetted:
	if(length(namesposf)>0){
		thisind=TRUE
		for(i in seq_along(namesposf)){
			thisind=thisind & do.call(ineq,list(min(range[[i]]),data[[namesposf[i]]])) & do.call(ineq,list(data[[namesposf[i]]],max(range[[i]])))
			}
		presentnames=intersect(c("psxf","psyf","pszf"),names(data))
		# Insert NA or subset the data:
		if(insert.NA){
			data[presentnames]=lapply(data[presentnames],function(x) {x[!thisind]=NA; x})
			}
		else{
			data[presentnames]=lapply(data[presentnames],function(x) x[thisind])
			}
		# Add to 'this$ind' (which are indices, not logical):
		this$ind[intersect(presentnames,treated)]=lapply(this$ind[intersect(presentnames,treated)],function(x) x[thisind])
		this$ind[setdiff(presentnames,treated)]=rep(list(which(thisind)),length(setdiff(presentnames,treated)))
		treated=union(treated,presentnames)
		}	
	
	# Acoustic data, simultaneously subsetted along with positions of the voxels:
	if(length(namesacoustic)>0){
		thisind=TRUE
		# Only one step in the for loop, but keeping the same structure as the above for loops for convenience:
		for(i in seq_along(namesacoustic)){
			thisind=thisind & do.call(ineq,list(min(range[[i]]),data[[namesacoustic[i]]])) & do.call(ineq,list(data[[namesacoustic[i]]],max(range[[i]])))
			}
		presentnames=intersect(c("psxx","psyx","pszx","vbsc","mvbs","volx","pr0s"),names(data))
		# Insert NA or subset the data:
		if(insert.NA){
			data[presentnames]=lapply(data[presentnames],function(x) {x[!thisind]=NA; x})
			}
		else{
			data[presentnames]=lapply(data[presentnames],function(x) x[thisind])
			}
		# Add to 'this$ind' (which are indices, not logical):
		this$ind[intersect(presentnames,treated)]=lapply(this$ind[intersect(presentnames,treated)],function(x) x[thisind])
		this$ind[setdiff(presentnames,treated)]=rep(list(which(thisind)),length(setdiff(presentnames,treated)))
		treated=union(treated,presentnames)
		}	
	
	# Other variables:
	for(i in seq_along(namesrange)){
		thisind=do.call(ineq,list(min(range[[i]]),data[[namesrange[i]]])) & do.call(ineq,list(data[[namesrange[i]]],max(range[[i]])))
		# Insert NA or subset the data:
		if(insert.NA){
			data[[namesrange[i]]][thisind]=NA
			}
		else{
			data[[namesrange[i]]]=data[[namesrange[i]]][thisind]
			}
		# Update 'this$ind' (which are indices, not logical):
		if(namesrange[i] %in% treated){
			this$ind[[namesrange[i]]]=this$ind[[namesposx[i]]][thisind]
			}
		else{
			this$ind[[namesrange[i]]]=which(thisind)
			}	
		}
	treated=union(treated,namesrange)
	
	
	##### Output #####
	# Add 'subsetind' is required:
	if(ind.out && length(this$ind)>0){
		names(this$ind)=paste("ind_",treated,sep="")
		data[["subs"]]=this$ind
		names(data[["subs"]])=treated
		}
	# Output:
	data[["treated"]]=treated
	data
}
