#*********************************************
#*********************************************
#' Extracts a subset of TSD data according to the array subset 'ind' and/or the cartesian subset 'range' and/or the logical/numeric vector of subscripts 'subset'. 
#'
#' @param data  is a list of elements named according to the TSD file format.
#' @param ind  is a list of indexes, as typed into the [] of an array, where 0 and NULL denotes all indexes.
#' @param range  is a list of elements with names matching names if 'data', specifying the range of the corresponding elements.
#' @param subset  is a numeric or logical vector/expression indicating elements or rows to keep (as used in subset()). Missing values are taken as false, and subset=0 or subset=NULL indicates no subsetting.
#' @param ind.out  is TRUE if the indexes for the elements segmented are to be returned.
#' @param drop  is TRUE if dimensions of only one level is to be removed from the output.
#' @param strict  is TRUE if strict inequality is to be used when subsetting according to 'range'.
#' @param insert.NA  is TRUE if the discarded data are to be kept as NA.
#' @param only.match  is TRUE if only the arrays of length equal to the length of 'subset' are to be subsetted using 'subset'.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD arr.ind2ind dim_all extractIndSubset ind.expand len_all NAs
#'
#' @export
#'
subset_TSD <- function(data, ind=list(), range=list(), subset=NULL, ind.out=FALSE, drop=TRUE, strict=TRUE, insert.NA=FALSE, only.match=FALSE, pad=c("end","start")){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2010-03-23 - Clean version.
	# Update: 2010-03-23 - Updated to only change the elements named "vbsc" or "mvbs, and "psx*", "psy*" and "psz*" (and oprionally "volx) Previously the other elements of the input 'data' were lost.
	# Update: 2010-05-10 - Method changed completely to a more robust and general, returning 'ind' if required. Also, the name of the array subscripts changed from 'subset' to 'ind' in accordance to the subset() function.
	# Update: 2011-01-17 - Added subsetting based on acoustic values.
	# Update: 2011-03-17 - Added the option 'affect' for specifying variables to subset using the result of 'range'.
	# Update: 2011-05-20 - Added the option 'only.match' for subsetting only the arrays of the same length as 'subset'.
	# Update: 2013-04-02 - Removed 'affect' and fixed bugs.
	# Last: 2013-08-07 - Implemented the function extract.range.TSD().
	########### DESCRIPTION: ###########
	# Extracts a subset of TSD data according to the array subset 'ind' and/or the cartesian subset 'range' and/or the logical/numeric vector of subscripts 'subset'. 
	########## DEPENDENCIES: ###########
	# extractIndSubset(), extract.range.TSD()
	############ VARIABLES: ############
	# ---data--- is a list of elements named according to the TSD file format.
	# ---ind--- is a list of indexes, as typed into the [] of an array, where 0 and NULL denotes all indexes.
	# ---range--- is a list of elements with names matching names if 'data', specifying the range of the corresponding elements.
	# ---subset--- is a numeric or logical vector/expression indicating elements or rows to keep (as used in subset()). Missing values are taken as false, and subset=0 or subset=NULL indicates no subsetting.
	# ---ind.out--- is TRUE if the indexes for the elements segmented are to be returned.
	# ---drop--- is TRUE if dimensions of only one level is to be removed from the output.
	# ---strict--- is TRUE if strict inequality is to be used when subsetting according to 'range'.
	# ---insert.NA--- is TRUE if the discarded data are to be kept as NA.
	# ---only.match--- is TRUE if only the arrays of length equal to the length of 'subset' are to be subsetted using 'subset'.
		
	
	##################################################
	##################################################
	##### Preparation #####
	# 'strict' is used when subsetting according to 'range':
	if(strict){
		ineq <- "<"
		}
	else{
		ineq <- "<="
		}
	
	# Input 'data' needs to be a list:
	if(!is.list(data)){
		stop("Input 'data' must be a list of elements named according to the TSD file format")
		}
	# Store the dimensions of the data:
	dims <- lapply(data, dim_all)
	# Insert lengths for lists:
	arelists <- sapply(data, is.list)
	if(any(arelists)){
		dims[arelists] <- len_all(dims[arelists])
		}
	# Store the number of dimensions of the data:
	ndims <- sapply(dims, length)
	
	# Input 'ind' needs to be a list:
	if(!is.list(ind)){
		ind <- list(ind)
		}
	# If the length of 'ind' is 1, a warning is issued:
	lind <- length(ind)
	#if(lind==1 && any(ndims>1)){
	#	warning("'ind' contains only one vector. All arrays of one or more dimensions will be subsetted along the first dimension")
	#	}
	# Do nothing if 'ind', 'range' and 'subset' are non-specified:
	treated <- names(data)[ndims>=lind]
	if(all(length(ind)==0, length(range)==0, length(subset)==0)){
		if(ind.out){
			data[["treated"]] <- treated
			data[["subs"]] <- lapply(data[treated], function(x) arr.ind2ind(ind.expand(list(),dim(x)),dim(x)))
			}
		return(data)
		}
	
		
	##### Execution and output #####
	# Extract the subsets specified by 'ind' and 'subset'. The data are affected, and the indices of the subsets are stored as well, to be intersected with subsets in the following:
	if(!any(identical(subset,0), length(subset)==0) || !any(identical(ind,0), length(ind)==0)){
		if(length(treated)==0){
			warning("'ind' has more elements than the dimensions of any of the elements of 'data'")
			this <- list()
			}
		if(length(treated)>0){
			# Get the maximum length in each dimension:
			maxdim <- dims[treated]
			nmaxdim <- max(ndims[treated])
			maxdim <- lapply(maxdim, function(xx) c(xx, NAs(nmaxdim-length(xx))))
			maxdim <- apply(as.data.frame(maxdim), 1, max, na.rm=TRUE)
			if(any(maxdim > sapply(ind, length)) || !any(identical(subset,0), length(subset)==0)){
				this <- extractIndSubset(data[treated], ind=ind, subset=subset, ind.out=TRUE, drop=drop, insert.NA=insert.NA, only.match=only.match, pad=pad)
				data[treated] <- this$x
				}
			else{
				this <- list()
				}
			}
		}
	else{
		this <- list()
		treated <- NULL
		}
	
	# Extract the subset given by 'range' and return:
	extract.range.TSD(data=data, range=range, this=this, treated=treated, ineq=ineq, ind.out=ind.out, insert.NA=insert.NA)
	##################################################
	##################################################
	}
