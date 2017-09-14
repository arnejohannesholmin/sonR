#*********************************************
#*********************************************
#' Extracts variables from acoustic events.
#'
#' @param event  is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.
#' @param var  is a string vector representing the variables to return. The function searches for the variables in the files contained in the event directory. The legal variable names are listed below:
#' @param t  is either a vector of the numbers of the pings to be returned, as listed from 1 to the number of pings in the event, or a vector of time points given as strings "yyyymmddHHMMSS.FFF" or "HHMMSS.FFF" from which the range of the time points to be read is extracted. If t == "all", all files are read and if t == "none" an empty list is returned.
#' @param cruise  is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event to be read.
#' @param ind  is a list of indexes, as typed into the [] of an array extracting a subset of the acoustic data (vbsc or mvbs).
#' @param TIME  is an optional list as returned from UNIX_time() (speed up reading).
#' @param adds  is an optional list of variables overriding the variables located in the event directory, used when calculating the positions and volumes of the voxels. Elements of 'adds' representing dynamic variables of the vessel (like "rtzv", "psxv" and so on) need to be at time steps 't' and thus have length equal to the length of 't'.
#' @param esnm  is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.
#' @param TVG  is FALSE if TVG compensation is to be removed from the data. No longer active. TVG.exp = 0 now corresponds to TVG = FALSE.
#' @param TVG.exp  is the exponent in the range compensation function (Time varied gain or TVG). TVG.exp=2 corresponds to the convensional 20 log R-TVG used for volume backscattering strength Sv. TVG.exp=0 removes the range compensation.
#' @param ideal  is TRUE if the speed of sound in water is assumed to be constant.
#' @param cs  is the coordinate system of the voxel midpoints or edges of the voxels.
#' @param seabed  is the depth of the seabed, at which the beams are reflected when calculating the midpoints of voxels.
#' @param rot  see soundbeam.TSD().
#' @param compensation  is a vector of string giving which rotation values that are compensated for in the sonar. Only c("pitch", "roll") is available for the current version. Used in soundbeam.TSD.
#' @param exact  is TRUE if the volumes should be calculated by the slower procedure incorporating the edges of the voxels. If FALSE (default) use the fast functions which only account for changes in the thickness of the voxels.
#' @param origin  is (1) a vector of two elements representing the origin of the global coordinate system (G), (2) the numbering index of the ping in the total sequence of pings of the event, which is to be regarded as the origin of (G) (ignoring heave so that the x-y-plane of (G) is on the surface of the sea), or (3) NULL, implying that the origin be put to the mid point of the vessel posistions.
#' @param dir.data  is the path to the directory in which the projects are stored, defaulted by the variable Acoustics_datasets_directory().
#' @param drop.out  is FALSE if output should not be dropped of dimensions of only one level.
#' @param strip.out  is FALSE if duplicated elements of the output list should be kept.
#' @param dir.out  is TRUE if only the path to the directory of the event is to be returned.
#' @param Paout  is TRUE if pressure data are to be returned in Pascal.
#' @param bgns  indicates wheter the estimated background noise should be subtracted, resulting in some negative sv-values. If TRUE (default) the background noise is kept in the data (no action). If FALSE the background noise estimates are attemptedly located in the event directory or the noise estimate directory of the cruise. If string, it should be the path to the file holding the background noise estimates. If the length of 'bgns' is longer than 1, the background noise estmate taken as the simple mean for all beams, regardless of the periodic noise is used.
#' @param pdns  indicates wheter the estimated periodic noise should be subtracted. If TRUE (default) the periodic noise is kept in the data (no action). If FALSE the periodic noise estimates are attemptedly located in the event directory or the noise estimate directory of the cruise. If string, it should be the path to the file holding the periodic noise estimates.
#' @param nrns  indicates wheter the estimated close range (previously named "near range") noise should be subtracted. If TRUE (default) the close range noise is kept in the data (no action). If FALSE the near range noise estimates are attemptedly located in the event directory or the noise estimate directory of the cruise.
#' @param hins  indicates wheter the estimated high intensity noise should be subtracted. If TRUE (default) the high intensity noise is kept in the data (no action). If FALSE the high intensity noise estimates are attemptedly located in the event directory or the noise estimate directory of the cruise.
#' @param kern  is the standard deviation in units of the number of voxels used in a Gaussian kernel smoother along the beams (no smoothing if kern == 0 or length(kern) == 0, which is the default). If given as an integer, say 5L, median smoothing is applied instead of Gaussian.
#' @param merge  is TRUE if variables read for more than one file are to be merged.
#' @param segpar  is a list of elements named "bwGp", "lsth"/"rlst", "usth"/"rust", or "sgth"/"sgt0" specifying the parameters of the segmentation data to read.
#' @param pamkpar  is a list of parameters ('krange', 'criterion', 'alpha' and 'mindist') used pamk() when clustering the segmented voxels 'sgsc'. The voxels 'sgsc' are also ordered so that the voxels belonging to the largest cluster lead and the second to largest cluster follows. A suggested set of parameters are pamkpar = list(krange = 1:4, criterion = "asw", alpha = 0.05, N = 1e2, mindist = 100).
#' @param nsind  is a vector of indexes along the beams, as input to ind.expand(), used to select the subset over which the estimation of the phase of the periodic noise is done. If given as a single numeric, the outermost 'nsind' voxels are used in each beam.
#' @param hins_add  is the number of voxels that should be discarded on both sides of high intensity noise voxels voxels along beams, used for accounting for possible high values that are related to the high intensity noise but not classified as such voxels.
#' @param pdns_scale  is used in get.pdns_phase.event() to scale the noise in order to allow the optimization to work.
#' @param TOV  is the time offset of the vessel information, discovered by Holmin and Korneliussen in 2013. The default is found in "/Applications/echoIBM/Documentation/Error in yaw MS70/Error in yaw MS70.R".
#' @param allow.old  is a TRUE if old UNIX_time file is accepted (still with the correct list of files).
#' @param onestep  is TRUE to allow for files with only one time step to be read regardless of 't'.
#' @param ...  are inputs used in ftim.TSD(), but more importantly in pplot3d.TSD(). Particularly, 'ind', 'range', and 'subset' can be used to subset the data in pplot3d.TSD(), but 'ind' can also be used to subset the acoustic data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom SimradRaw apply.TVG
#' @importFrom TSD arr.ind2ind dim_all ftim.TSD ftim2utim ind.expand info.TSD labl.TSD mergeListKeepDimensions mtim.TSD ones read.TSD read.TSDs strff utim2mtim zeros prettyIntegers
#' @importFrom tools file_ext
#' @importFrom data.table rbindlist
#' @importFrom stats setNames
#'
#' @export
#' @rdname read.event
#'
read.event<-function(event=1, var="pings", t=1, cruise=2009116, TIME=FALSE, adds=NULL, esnm="MS70", TVG=TRUE, TVG.exp=2, ideal=TRUE, cs="g", seabed=-12000, rot=2, compensation=c("pitch", "roll"), exact=FALSE, origin=1, dir.data=NULL, drop.out=TRUE, strip.out=TRUE, dir.out=FALSE, mask=TRUE, Paout=TRUE, bgns=TRUE, pdns=TRUE, nrns=TRUE, hins=TRUE, kern=NULL, merge=TRUE, msg=TRUE, segpar=NULL, pamkpar=list(), nsind=0.75, hins_add=10, pdns_scale=1e-14, TOV=0, allow.old=FALSE, ggsz=10, pad=TRUE, split=TRUE, cal=1, fanWidth="b2", onestep=TRUE, other=FALSE, fill=NA, ...){
	
	############ AUTHOR(S): ############
	# Arne Johannes Holmin
	############ LANGUAGE: #############
	# English
	############### LOG: ###############
	# Start: 2009-03-10 - First version.
	# Update: 2009-06-04 - Cleaned up and documented.
	# Update: 2009-09-04 - Changed corresponding to the change in HIB(), where integer() as output represents no HIBs found.
	# Update: 2010-02-11 - Method changed completely to support "cases" and the TSD file format.
	# Update: 2010-03-10 - Added support for specifying time range by 't', given as yyyymmddHHMMSS.FFF pr HHMMSS.FFF.
	# Update: 2010-03-22 - Changed for loop through time steps for reading dynamic school variables to beeing av for loop through files (speeding up function).
	# Update: 2010-04-17 - Added support for reading rtxf and rtzf when only velocity data are present.
	# Update: 2010-05-04 - Added option 'pingind' for returning the list 'p' used when reading sv data.
	# Update: 2010-05-17 - Changed to read school data.
	# Update: 2010-11-15 - Added the option 'indt'.
	# Update: 2010-11-17 - Fixed bug when var = "school".
	# Update: 2010-11-17 - Fixed the reading when var = "school" to using read.TSDs().
	# Update: 2011-01-08 - Major changes: Mehtod of reading altered and simplified to first establishing a list of time indexes for each file, and then use this when reading the data. Also, removed the option 'indt', requiring that time steps given by 't' are always general time step indexes.
	# Update: 2011-01-08 - Removed extra merging of the school data which was done because of a misunderstanding of the method in read.TSDs(). Merging is done in read.TSDs().
	# Update: 2011-02-19 - Added the option 'ts.show', used for displaying information about the time steps of the files located in the event. Also changed the default value og 'merge' to TRUE.
	# Update: 2011-03-08 - Significant update: When reading time variables ("time", "mtim", "utim", "ftim", "indt") these variables are only read from files that are not ctd-filea or beams-files, because the time information in these files is not relevant. Also the time variables are sorted and made unique.
	# Update: 2011-03-09 - Improvement of the time used when reading events of many time steps. The option of reading the file "UNIX_time.tsd" is added, and writing the file if not present. This file contains the unix time points and file names (without path) of the files present in the event, and reading this is much faster than reading unix time for all files (done by a for loop as long as the number of time steps).
	# Update: 2011-05-19 - Fixed a bug when static variables are present in the dynamic school files. The former variable 'schooltype' is replaced by the two variables 'schooltypeD' and 'schooltypeS', so that a file can be read as both dynamic and static school file.
	# Update: 2011-05-20 - Added the option 'strip.out'.
	# Update: 2011-06-28 - Changed to detect ctd and school files using variables, and not only file estensions (in order to read beam pattern files with the ".tsd" file extension).
	# Update: 2011-08-08 - Changed the treatment of the time information to discard the time points stored in .beams files and .ctd files.
	# Update: 2011-08-11 - Added the function UNIX_time to the treatment of time. This involves writing new UNIX_time files for all that miss the 'i000' variable.
	# Update: 2011-09-08 - Added the option 'ind' for extracting subsets of acoustic data.
	# Update: 2011-09-19 - Changed the method of removing HIBs in the acoustic data to using the function rm.HIB().
	# Update: 2011-09-25 - Changed to using volx.TSD() when calculating voxel volumes, and reading beams-data, ctd-data and vessel-data if any voxel data are requested.
	# Update: 2011-10-20 - Added the option 'cs' used in soundbeam.TSD() for generating voxel position values.
	# Update: 2011-10-22 - Changed the effect of 'merge' for time variables so that if all time points are present in the vessel file, these are returned, else all time variables of the same type ('mtim', 'indt' ...) are merged.
	# Update: 2011-10-22 - Added the function getSchoolfileType().
	# Update: 2012-02-20 - Changed to use the function echoIBM.getSchoolfileType().
	# Update: 2012-07-01 - Added the function event.path() for extracting the path to the event from the parameters 'event' and 'cruise'.
	# Update: 2012-07-25 - Added 'segpar' to the parameter list, for selecting segmentation data by specifying the segmentation parameters used.
	# Update: 2012-07-25 - Radical change in the way "UNIX_time.TSD" is written, and the information used in the function to reduce CPU time. It was observed that the files are read redundently to locate the information that is already present in "UNIX_time.TSD". The change in UNIX_time() involves discarding beams-files and ctd-files when determining the time indexes 'i000'. The changes in read.event() involves (1) not reading 'indt' for all files, which was done previously to insert 'indt' directly from the files which included this information, (2) using 'r000' as extracted from "UNIX_time.TSD" in read.TSD() and read.TSDs(), (3) numerous other changes which speeds up the function radically.
 	# Update: 2012-08-01 - Another radical change. Separated much of the function into sub-functions. Restructured to reduce CPU time.
 	# Update: 2012-08-01 - Changed 'tlist' as input to 'TIME' as input, and similarly for the optional output.
 	# Update: 2012-08-06 - Fixed bug when reading segmentation data.
	# Update: 2012-10-12 - Changed to use the function list.files_caseInsensitive() instead of list.files(), due to the difference in the performance of the latter function in the terminal version of R and in the GUI version. The custom function list.files_caseInsensitive() sorts in a case insensitive matter, as is done by the GUI version of R.
	# Update: 2013-01-05 - Removed the subtraction of high intensity beams, as the high intensity noise is accounted for by the method presented in the second paper of the PhD of Holmin.
	# Update: 2013-03-21 - Removed info in UNIX_time().
	# Update: 2013-05-06 - Added using correctVessel() to correct for time offset of the vessel information.
	# Update: 2013-05-07 - Added the incertion of NAs for variables that are only present in a subset of the time steps, and which are not of type "beams", "ctd", "vbsc" and the like, "vessel", "voxels", and "time".
	# Update: 2013-08-09 - Added origin == NULL, implying that the origin be put to the mid point of the vessel posistions. Also, simplified the variable list.
	# Update: 2013-10-01 - Added the function removeDuplicated_tlist() which merges the time list to discard duplicated time steps.
	# Update: 2013-10-09 - Included removeDuplicated_tlist() in read.event_read_sgsc() in order to avoid warnings associated with duplicated time steps before the selection of segmentation files is done.
	# Update: 2014-05-05 - Added reading 'vxIX' to decompress the acoustic data.
	# Update: 2014-09-30 - Removed 'r000' accoring to the changes in read.TSD() and write.TSD().
	# Update: 2014-09-30 - Fixed bugs in reading "vbsC" and "vbsA".
	# Update: 2015-04-24 - Changed to use for loop when calculating the voxel positions psxx, psyx, and pszx.
	# Last: 2015-11-02 - Disactivating 'TVG', and only using 'TVG.exp'.
	########### DESCRIPTION: ###########
	# Extracts variables from acoustic events.
	########## DEPENDENCIES: ###########
	# read.TSD(), global2car(), soundbeam.TSD(), volx.MS70(), HIB(), apply.TVG(), zeros(), utim.TSD(), read.TSDs(), UNIX_time(), labl.TSD()
	############ DETAILS: ############
	# info.TSD() attempts to locate a table named "TSD_description_table" with two columns, where the first column contains the TSD variable names, and the second the descriptions of the variables. If not present, read.TSD_description_table() tries to read the table from the file given by 'file'. The used should create a file with two columns separated by "\t" with the variable names and desciptions as described above, with a header line "Label\tDescription". To avoid that this file is read each time info.TSD() is run, run the following code at the start of the R-session, where 'file' is the full path to the file of the TSD names and descriptions:
	# TSD_description_table = read.TSD_description_table(file)
	############ VALUE: ############
	# 
	############ REFERENCES: ############
	#
	############ SEAALSO: ############
	#
	############ EXAMPLES: ############
	# # Create a character matrix 'TSD_description_table' if not alreaddy existing:
	# if(!exists("TSD_description_table")) TSD_description_table = cbind(c("var1","var2"), c("First variable","Second variable"))
	# x <- list(var1 = 1:12)
	# info.TSD(x)
	############ VARIABLES: ############
	#
	#	Variable identifiers reading whole files or groups of variables (see the "TSD-format_names.dat" for all currently documented TSD variables):
	#		"pings" reads the variables located in .pings files
	#		"beams" all beams variables: labl.TSD("b")
	#		"beams0" reads (mostly) relevant beams variables: labl.TSD("rb")
	#		"vessel" reads vessel variables: labl.TSD("v")
	#		"ctd" CTD (conductivity, temperature, depth) data: labl.TSD("ctd")
	#		"voxels" reads voxel position and volume data: labl.TSD("vx")
	#		"time" reads time variables: labl.TSD("t")
	#		"school" read static and dynamic school variables: labl.TSD(c("ss", "ds"))
	
	
	##################################################
	##################################################
	##### Preparation #####
	filetypes = c("school", "beams", "ctd", "vessel", "pings", "seg", "tsd")
	filetypeslist = setNames(as.list(filetypes), filetypes)
	
	##########
	# Store the arguments to be passed on to pplot3d.TSD():
	ll = list(...)
	# Functions for reading noise data:
	getbgns = function(bgns){
		if(length(bgns)>1){
			type = 0
			bgns = bgns[1]
		}
		else{
			type = "s"
		}
		varname = paste("bgn", type, sep="")
		suppressWarnings(read.TSDs(bgns, var=varname, msg=FALSE))
	}
	getpdns = function(pdns, relevantpdnsvar){
		suppressWarnings(read.TSDs(pdns, var=relevantpdnsvar, msg=FALSE))
	}
	getnrns = function(nrns, mode=c("p", "a")){
		varname = paste("nr0", substr(mode[1], 1, 1), sep="")
		suppressWarnings(read.TSDs(nrns, var=varname, msg=FALSE))
	}
	gethins = function(hins){
		suppressWarnings(read.TSDs(hins, var=labl.TSD("hins"), msg=FALSE))
	}
	if(isTRUE(TIME)){
		TIMEout = TRUE
	}
	else{
		TIMEout = FALSE
	}
	##########
	
	
	##########
	# If not given, set the data directory containing data (structured in the directory structure specified in the documentation of echoIBM) as the string Acoustics_datasets_directory():
	if(is.null(dir.data)){
		dir.data = Acoustics_datasets_directory()
	}
	# Locate the event:
	if(length(event) && identical(file.info(event[1])$isdir, FALSE)){
		filelist = event
		event=dirname(event[1])
	}
	else{
		filelist = NULL
	}
	event = event.path(event=event[1], cruise=cruise, esnm=esnm, dir.data=dir.data)
	esnm = event$esnm
	event = event$event
	
	# If the event does not exist, NULL is returned:
	if(!file.exists(event)){
		warning(paste("Event ", event, " missing", sep=""))
		return(list())
	}
	# If dir.out == TRUE, only the string of the tsd directory is returned (mainly for use in HIB.TSD()):
	if(dir.out){
		return(event)
	}
	##########
	
	
	##########
	# Accept strings with leading "-" or "+":
	pm = substr(var, 1, 1) %in% c("+", "-")
	if(any(pm)){
		var[pm] = substr(var[pm], 2, 5)
	}
	
	# If regenerated uniformly distributed points reflecting the vbsc, or if the volume backscattering coefficient in grid cells ("vbsC") are requested, read also voxel and vessel variables:
	psr = c("psxr", "psyr", "pszr")
	readpsr = any(var %in% c(psr, "vbsC", "vbsA"))
	if(readpsr){
		var = c("vbsc", "vessel", "voxels", "beams0", var)
	}
	
	# Are sv values to be read:
	readsv = any(var %in% c("vbsc", "mvbs", "pings", psr))
	# Define various variable chategories:
	staticschoolnames = labl.TSD("ss")
	dynschoolnames = labl.TSD("ds")
	beamsnames = labl.TSD("b")
	ctdnames = labl.TSD("ctd")
	vesselnames = labl.TSD("v")
	relevantbeamsvar = labl.TSD("rb")
	relevantctdvar = labl.TSD("rc")
	relevantpdnsvar = labl.TSD("rp")
	segvar = labl.TSD("sg")
	voxelvar = labl.TSD("vx")
	#timevar = labl.TSD("t")
	timevar = c("utim","mtim","indt")
	
	# Insert school variables and voxel variables to 'var', and store the unaltered 'var':
	oldvar = var
	if("school" %in% var){
		var = setdiff(c(var, staticschoolnames, dynschoolnames), "school")
	}
	if("voxels" %in% var){
		var = setdiff(c(var, voxelvar), "voxels")
	}
	if("beams0" %in% var){
		var = setdiff(c(var, relevantbeamsvar), "beams0")
	}
	if("seg" %in% var){
		var = setdiff(c(var, segvar), "seg")
	}
		if("time" %in% var || TOV != 0){
		timevar_present = timevar
	}
	else{
		timevar_present = intersect(timevar, var)
	}
	
	# Only the header is read if t == "none" or var == "none" (or none of the elements of 'var' are available):
	if(identical(t, "none") || identical(var, "none")){
		var = NULL
		t = 0
	}
	if(length(t) == 0 || is.na(t)){
		return()
	}
	##########
	
	
	##########
	# The list of files in the tsd-directory inside the case:
	if(length(filelist) == 0){
		# Using recursive=TRUE discards empty directories:
		filelist = list.files_caseInsensitive(event, full.names=TRUE, recursive=TRUE)
	}
	if(length(filelist) == 0){
		warning("The event is empty")
		return(list())
	}
	
	# Remove the file "UNIX_time.tsd" from the list if it is present, as this only contains information about the unix time points of the other files:
	#filelist = filelist[basename(filelist) != "UNIX_time.tsd"]
	filelist = rmUnix_Sgpm(filelist)
	# Selecting only files with valid extensions:
	ext = tolower(file_ext(filelist))
	fInd = lapply(filetypeslist, function(xx) which(ext==xx))
	
	# Discard the files that are not relevant to the requested variables:
	relevantFiles = NULL
	# 'tempvar' is only used here to reduce the number of files read.
	tempvar = var
	
	if("seg" %in% var || any(segvar %in% var)){
		relevantFiles = c(relevantFiles, fInd$seg)
		tempvar = setdiff(tempvar, segvar)
	}
	# If any time variables are requested, all files are relevant:
	if(any(c(timevar, "time") %in% timevar_present)){
		# If a clean list of time steps is requested, only read one time step, which will include the 'I000' and 'U000' variables:
		if(merge){
			relevantFiles = c(relevantFiles, fInd$vessel)
		}
		else{
			relevantFiles = c(relevantFiles, seq_along(filelist))
		}
		tempvar = setdiff(tempvar, c(timevar, "time"))
	}
	if("beams" %in% var || any(var %in% beamsnames) && length(fInd$beams)>0){
		relevantFiles = c(relevantFiles, fInd$beams)
		tempvar = setdiff(tempvar, c("beams", beamsnames))
	}
	if("ctd" %in% var || any(var %in% ctdnames) && length(fInd$ctd)>0){
		relevantFiles = c(relevantFiles, fInd$ctd)
		tempvar = setdiff(tempvar, c("ctd", ctdnames))
	}
	# Allways keep the vesselfiles in order to read the vessel information (needed if 't' is given as "HHMMSS" and when voxel positions are needed)
	relevantFiles = c(relevantFiles, fInd$vessel)
	tempvar = setdiff(tempvar, "vessel")
	if(any(var %in% vesselnames) && length(fInd$vessel)>0){
		relevantFiles = c(relevantFiles, fInd$vessel)
		tempvar = setdiff(tempvar, vesselnames)
	}
	if(any(c("psxx", "psyx", "pszx") %in% var) || any(c("volx", "harx") %in% var)){
		relevantFiles = c(relevantFiles, fInd$beams, fInd$ctd, fInd$vessel)
		tempvar = setdiff(tempvar, voxelvar)
	}
	if("seg" %in% var || any(segvar %in% var)){
		relevantFiles = c(relevantFiles, fInd$seg)
		tempvar = setdiff(tempvar, segvar)
	}
	if(readsv){
		relevantFiles = c(relevantFiles, fInd$pings, fInd$beams, fInd$ctd)
		tempvar = setdiff(tempvar, c("vbsc", "mvbs", "pings", "psxr", "psyr", "pszr"))
	}
	if(any(c("angl", "angt") %in% var)){
		relevantFiles = c(relevantFiles, fInd$pings)
		tempvar = setdiff(tempvar, c("angl", "angt"))
	}
	if(any(var %in% c(staticschoolnames, dynschoolnames))){
		relevantFiles = c(relevantFiles, fInd$school)
		tempvar = setdiff(tempvar, c(staticschoolnames, dynschoolnames))
	}
	if(length(tempvar)>0){
		relevantFiles = c(relevantFiles, c(fInd$vessel, fInd$beams, fInd$ctd, fInd$school, fInd$seg, fInd$tsd))
	}
		
	# Add the files with file extension tsd, if TOV != 0:
	if(TOV != 0 && strff("ms70", esnm)){
		relevantFiles = c(relevantFiles, fInd$tsd)
	}
	
	# Discard the irrelevant files, and re-extract the file extensions and file types:
	relevantFiles = unique(relevantFiles)
	# Selecting only files with valid extensions:
	filelist = filelist[relevantFiles]
	ext = tolower(file_ext(filelist))
	fInd = lapply(filetypeslist, function(xx) which(ext==xx))
	##########
	
	
	##########
	# Read UNIX_time file if not given as input:
	if(is.list(TIME)){
		TIME[c("f000", "i000", "u000", "l000", "t000")] = lapply(TIME[c("f000", "i000", "u000", "l000", "t000")],function(x) x[relevantFiles])
	}
	else{
		TIME <- UNIX_time(event=event, file=TRUE, var="all", t=relevantFiles, msg=msg, allow.old=allow.old)
		# Sone check for the success of reading the UNIX time (unclear why):
		if(length(TIME$I000[[1]]) == 0 || length(TIME$U000[[1]]) == 0){
			TIME <- UNIX_time(event=event, file=TRUE, fresh=TRUE, msg=msg, allow.old=allow.old)
		}
	}
	# Rename to allow ftim.TSD() and mtim.TSD():
	if(!"utim" %in% names(TIME)){
		names(TIME)[names(TIME) == "u000"]="utim"
	}
	if(!"indt" %in% names(TIME)){
		names(TIME)[names(TIME) == "i000"]="indt"
	}
	if(!"labl" %in% names(TIME)){
		names(TIME)[names(TIME) == "l000"]="labl"
	}
	# Repport an error if there are no time steps in the event:
	if(sum(unlist(lapply(TIME$indt,length))) == 0){
		warning("The event \"", event, "\" has no time steps")
		return(list())
	}
	##########
	
	
	##########
	# Get the number of unique time point:
	numt = length(TIME$I000)
	
	# 't' may be given as formated time "yyyymmddHHMMSS.FFF" or "yyyymmddSSSSS.FFF":
	nchart = nchar(t)
	if(is.character(t) && length(nchart)>0 && min(nchart)>5){
		# Extract formated time if 't' is given as formated time:
		if(!"ftim" %in% names(TIME) && "ftim" %in% timevar_present){
			TIME$ftim = ftim.TSD(TIME[setdiff(names(TIME),"indt")])
		}
			
		t = as.character(t)
		# If the time is given in the classic as.POSIXlt format "yyyy-mm-dd HH:MM:SS", or with any other separating characters ":", "/", "-" or " ", strip 't' of these, and continue:
		if(length(grep(":", t))>0){
			t = strsplit(t, ".", fixed=TRUE)
			for(i in seq_along(t)){
				# Strip 't' of the separating characters, assuming that the time is given in the order of year, month, day, hour, minute, second:
				t[[i]][1] = gsub("[[:punct:] || [:blank:]]", "", t[[i]][1])
				#t[[i]][1] = gsub(":", "", t[[i]][1])
				#t[[i]][1] = gsub(";", "", t[[i]][1])
				#t[[i]][1] = gsub("/", "", t[[i]][1])
				#t[[i]][1] = gsub("-", "", t[[i]][1])
				#t[[i]][1] = gsub(" ", "", t[[i]][1])
				t[[i]] = paste(t[[i]], collapse=".")
			}
			t = unlist(t)	
		}
		# Constructing the unix time range based on 't':
		nchart = nchar(strsplit(t, ".", fixed=TRUE)[[1]][1])
		if(length(nchart)>0 && all(nchart == 6)){
			day = substr(read.TSD(filelist[fInd$vessel[1]], var="ftim", t=1)$ftim, 1, 8)[[1]]
			userutim = ftim2utim(paste(day, t, sep=""))
		}
		else if(length(nchart)>0 && all(nchart == 14)){
			userutim = ftim2utim(t)
		}
		else{
			stop("Invalid input 't'. Must be integers, \"none\", \"all\" or a 6-digit or 12 digit numbers giving the time as \"HHMMSS.FFF\" or \"yyyymmddHHMMSS.FFF\"")
		}
		# Selecting only the time points within the range of 't', or the closest if only one time point is given:
		if(length(userutim) == 1){
			t = which.min(abs(TIME$U000-userutim))
		}
		else{
			t = which(TIME$U000>min(userutim) & TIME$U000<max(userutim))
		}
		cat("Time points nr.", prettyIntegers(t), "selected\n")
	}
	else if(!is.numeric(t) && !identical(t, "all")){
		warning("Time input 't' must be numeric or \"none\", \"all\" or as formated time \"yyyymmddHHMMSS.FFF\" or \"yyyymmddSSSSS.FFF\" ('t' set to 0)")
		var = 0
	}
	# If t == "all", all time points are read:
	if(identical(t, "all")){
		t = TIME$I000
	}
	# Give a warning and discard any non-positive time step indices:
	if(any(t<= 0)){
		t = t[t>0]
	}
	
	# 'tlist' is a list of one vector for each file in 'filelist' holding the time step number in the files:
	tlist = vector("list", length(filelist))
	for(i in seq_along(filelist)){
		thisindt = match(as.numeric(t), TIME$indt[[i]])
		thisindt = thisindt[!is.na(thisindt)]
		if(length(thisindt)>0){
			tlist[[i]] = thisindt[!is.na(thisindt)]
		}
	}
	# Add ones to 'tlist' for the files with only one time step:
	if(onestep){
		onestep = unlist(lapply(TIME$indt, function(x) length(x) == 1)) & ext != "pings"
	}
	tlist[onestep] = as.list(ones(sum(onestep)))
	
	# Supply the 'fInd$school' with files that contain relevant variables, as these may be contained in other files than the ones with ".school" as file extension:
	if(!all(onestep) && any(var %in% c(staticschoolnames, dynschoolnames))){
		# For loop through the files:
		for(i in seq_along(filelist)){
			if(any(c(staticschoolnames, dynschoolnames) %in% TIME$labl[[i]])){
				fInd$school = c(fInd$school, i)
			}
		}
		# Clean the schoolfile indexes:
		fInd$school = unique(fInd$school)
	}
	
	# Keep only files with positive length in 'tlist', and re-extract the file extensions and file types:
	relevantFiles = which(sapply(tlist, length)>0)
	# Selecting only files with valid extensions:
	tlist = tlist[relevantFiles]
	filelist = filelist[relevantFiles]
	ext = tolower(file_ext(filelist))
	fInd = lapply(filetypeslist, function(xx) which(ext==xx))
	
	
	# Merging the list of .vessel-files and other files to obtain the list of files that are not .pings-files so that the .vessel-files are sorted first in the list. (At the end of the function list elements of duplicated names are removed so that only the first of the duplicated elements is chosen.)??? Why are vessel files mentioned?
	fInd$notpings = unlist(fInd[names(fInd)!="pings"])
	#fInd$notpings = unique(match(fInd$notpings, filelist))
	
	
	##### Execution and output #####
	out = list()
	#############################################
	### (0) Reading the number of time steps: ###
	#############################################
	# The number of time steps is alreaddy read:
	if("numt" %in% var){
		out$numt = numt
		# Remove "numt" from 'var':
		var = setdiff(var, "numt")
	}
	#############################################
	#############################################
	
	
	###################################
	### (1) Reading time-variables: ###
	###################################
	# Read time variables, which will be sorted and unique if merge == TRUE. Also ctd-files and beams-files will not be read for time (as often these files are used along with other files in a way that does not make the time of the ctd-files and the beams-files relevant):
	### (1a) Read indt:
	if("indt" %in% timevar_present){
		if(merge){
			out$indt = TIME$I000[t[t <= length(TIME$U000)]]
		}
		else{
			out$indt = TIME$indt
		}
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("time", "indt"))
	}
	### (1b) Read utim:	
	if("utim" %in% timevar_present){
		if(merge){
			out$utim = TIME$U000[t[t <= length(TIME$U000)]]
		}
		else{
			out$utim = TIME$utim
		}
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("time", "utim"))
	}
	### (1c) Read mtim:
	if("mtim" %in% timevar_present){
		if(merge){
			out$mtim = utim2mtim(TIME$U000[t[t <= length(TIME$U000)]])
		}
		else{
			out$mtim = mtim.TSD(TIME)
		}
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("time", "mtim"))
	}
	### (1d) Read ftim:
	if("ftim" %in% var){
		if(merge){
			#out$ftim = utim2ftim(TIME$U000[t])
			out$ftim = ftim.TSD(list(utim=TIME$U000[t[t <= length(TIME$U000)]]), ...)
		}
		else{
			out$ftim = ftim.TSD(TIME, ...)
		}
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("time", "ftim"))
	}
	###################################
	###################################

	####################################
	### (2) Reading beams-variables: ###
	####################################
	### (2a) Read all the beams-files:
	beamsvar = NULL
	if("beams" %in% var){
		beamsvar = "all"
		var = setdiff(var, "beams")
	}
	else if(any(var %in% beamsnames) && length(fInd$beams)>0){
		beamsvar = intersect(var, beamsnames)
	}
	if(length(beamsvar)){
		temp <- read.event_read_files(files=filelist, filesind=fInd$beams, tlist=tlist, var=beamsvar)
		#out[names(temp)] = temp
		addnames <- setdiff(names(temp), names(out))
		out[addnames] = temp[addnames]
		var = setdiff(var, beamsvar)
	}
	####################################
	####################################
	
	##################################
	### (3) Reading ctd-variables: ###
	##################################
	### (3a) Read all the ctd-files:
	if("ctd" %in% var){
		suppressWarnings(temp <- read.TSDs(filelist[fInd$ctd], msg=FALSE))
		#out[names(temp)] = temp
		addnames <- setdiff(names(temp), names(out))
		out[addnames] = temp[addnames]
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, "ctd")
	}
	### (3b) Read other ctd-variables:
	if(any(var %in% ctdnames) && length(fInd$ctd)>0){
		temp = read.event_read_other_ctd(var, ctdfiles=filelist[fInd$ctd], ctdfilesind=fInd$beams, TIME=TIME)
		#out[names(temp)] = temp
		addnames <- setdiff(names(temp), names(out))
		out[addnames] = temp[addnames]
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, ctdnames)
	}
	##################################
	##################################
	
	#####################################
	### (4) Reading vessel-variables: ###
	#####################################
	### (4a) Read all vessel information:
	vesselvar = NULL
	if("vessel" %in% var){
		vesselvar = "all"
		var = setdiff(var, "vessel")
	}
	else if(any(var %in% vesselnames) && length(fInd$vessel)>0){
		vesselvar = intersect(var, vesselnames)
	}
	if(length(vesselvar)){
		temp = read.event_read_files(files=filelist, filesind=fInd$vessel, tlist=tlist, var=vesselvar, origin=origin)
		#out[names(temp)] = temp
		addnames <- setdiff(names(temp), names(out))
		out[addnames] = temp[addnames]
		var = setdiff(var, vesselvar)
	}
	
	# If the time offset 'TOV' is given different from 0, extract any of "rtzv", "lonv", "latv", or "ispv" located in the raw dynamic vessel file:
	if(TOV != 0 && strff("ms70", esnm)){
		correctout = correctVessel(names(out), out[c("mtim", "utim")], filelist[fInd$tsd], TOV=TOV)
		out[names(correctout)] = correctout
	}
	#####################################
	#####################################
	
		
	################################################
	### (5) Reading voxel positions and volumes: ###
	################################################
	### (5a) Get the required vessel information:
	# Read vessel dynamic variables, beams variables and ctd-data, in the case that voxels positions (psxx, psyx, pszx) or voxel volumes (volx) are requested:
	if(any(c("psxx", "psyx", "pszx", "rngx", "volx", "harx") %in% var)){
		if("vessel" %in% oldvar){
			vessel = out[vesselnames]
		}
		else{
			vessel = read.event_read_files(files=filelist, filesind=fInd$vessel, tlist=tlist, var="all", origin=origin)
		}
		# If the time offset 'TOV' is given different from 0, extract any of "rtzv", "lonv", "latv", or "ispv" located in the raw dynamic vessel file:
		if(TOV != 0 && strff("ms70", esnm)){
			correctout = correctVessel(names(vessel), out[c("mtim", "utim")], filelist[fInd$tsd], TOV=TOV)
			vessel[names(correctout)] = correctout
		}
		
		# Read beams-files and ctd-files:
		if(all(relevantbeamsvar %in% names(out))){
			beams = out[relevantbeamsvar]
		}
		else{
			beams = read.event_read_files(files=filelist, filesind=fInd$beams, tlist=tlist, var="all")
		}
		if(all(relevantctdvar %in% names(out))){
			ctd = out[relevantctdvar]
		}
		else{
			ctd = read.TSDs(filelist[fInd$ctd], var=relevantctdvar, t="all", msg=FALSE)
		}
	}
	### (5b) Read the voxels positions:	
	if(any(c("psxx", "psyx", "pszx") %in% var)){
		# Midpoints of voxels:
		out[c("psxx", "psyx", "pszx")] <- psx.TSD(c(adds, vessel, beams, ctd), cs=cs, seabed=seabed, rot=rot, compensation=compensation, ideal=ideal, drop.out=drop.out, pad=pad, split=split)
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("psxx", "psyx", "pszx"))
}
	### (5c) Read the voxel volumes:
	if(any(c("volx", "harx") %in% var)){
		volxharxrequested = intersect(c("volx", "harx"), var)
		out[volxharxrequested] = volx.TSD(c(adds, beams, ctd, vessel), esnm=esnm, var=var, fanWidth=fanWidth)[volxharxrequested]
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("volx", "harx"))
}
	### (5c) Read the ranges to the voxel centers (only for one beam):
	if("rngx" %in% var){
		out$rngx <- soundbeam_range(beams, pos="mid")
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, "rngx")
}
	################################################
	################################################
	
	
	######################################
	### (6) Reading segmentation data: ###
	######################################
	# Take special care of the segmentation data, which may exist for different parameter specifications:
	if(any(segvar %in% var)){
		# Special case where only 'sgPM' is requested:
		if(sum(segvar %in% var) == 1 && "sgPM" %in% var){
			out["sgPM"] = list(sgPM(filelist[fInd$seg]))
	}
		else{
			# Read the segementation data:
			thisout = read.event_read_sgsc(filelist=filelist, var=c(var,"indt"), filesind=fInd$seg, segpar=segpar, TIME=TIME, tlist=tlist, merge=merge, msg=msg)
			# Add empty list elements where no data were present:
			if(merge){
				thisout = read.event_insertNA(thisout, TIME, fInd$seg[thisout$segfilenr], filelist, t, tlist, var, segvar)
			}
			out[names(thisout)] = thisout
			
			# Get clusters:
			if("clsz" %in% var){
				thisout = read.event_get.clsz(out["sgsc"], filelist[fInd$beams], fInd$beams, tlist, pamkpar, esnm)
				out$sgsc = thisout$sgsc
				out$clsz = thisout$clsz
			}
		}
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, segvar)
	}
	######################################
	######################################
			
	
	##################################
	### (7) Reading acoustic data: ###
	##################################
	# Read the relevant pings files to obtain acoustic data. The method in the following block may be changed to a simpler method in the future, inspired by how the school-files are read below:
	if(readsv){
		# Only read if 't' is not empty:
		if(length(t)>0){
			# We need the number of beams 'numb':
			#suppressWarnings(numb <- read.TSD(filelist[fInd$pings][1], var=c("numb", "lenb"), t=1, dimension=FALSE))
			#lenb = numb$lenb
			#numb = numb$numb
			
			# Read the acoustic data, and the voxel indices "vxIX" giving the voxels holding data if present:
			suppressWarnings(temp <- read.TSDs(filelist[fInd$pings], t=tlist[fInd$pings], var=c("vbsc", "vxIX", "numb", "lenb", "utim"), dimension=TRUE, merge=merge, indt=FALSE, drop.out=FALSE, msg=FALSE))
			
			#if(length(temp$vbsc)>0){
			#	out$vbsc = temp$vbsc
			#}
			#if(length(temp$vxIX)>0){
			#	out$vxIX = temp$vxIX
			#}
			# Unzip compressed vbsc:
			tempNotPresentInOut <- !names(temp) %in% names(out)
			out[names(temp)[tempNotPresentInOut]] <- temp[tempNotPresentInOut]
			out <- read.event_unzip_vbsc(out, pad=pad, split=split, t=t, fill=fill)
			
			
			# # If the voxel indices of compressed acosutic data are present, uncompress:
			# if(length(temp$vxIX)>0){
			# 	out <- read.event_unzip_vbsc(out, pad=pad, split=split, t=t, fill=fill)
			# }
			# # Otherwise simply return the acsoutic data:
			# else{
			# 	
			# }
			# 
			
			### # Filling in for compressed acoustic data:
			### if(length(out$vbsc)>0 && length(out$vxIX)>0){
			### 	# In the unlikely event that all time steps have equal number of positive voxels in the compressed mode, split into a list:
			### 	if(!is.list(out$vbsc)){
			### 		thisnumt = ncol(out$vbsc)
			### 		#out$vbsc = split(out$vbsc, rep(seq_len(thisnumt), each=length(thisnumt)/thisnumt))
			### 		#out$vxIX = split(out$vxIX, rep(seq_len(thisnumt), each=length(thisnumt)/thisnumt))
			### 		out$vbsc = as.data.frame(out$vbsc)
			### 		out$vxIX = as.data.frame(out$vxIX)
			### 	}
			### 	# Create a list of vbsc, where each ping is represented in each element of the list, and fill inn the non-empty values indicated by 'vxIX':
			### 	for(i in seq_along(out$vbsc)){
			### 		tempvbsc = NAs(max(temp$lenb[,i]), temp$numb[i])
			### 		tempvbsc[out$vxIX[[i]]] = out$vbsc[[i]]
			### 		out$vbsc[[i]] = tempvbsc
			### 	}
			### }
			### 	
			### # Rearranging sv-values in a 3 dimensional array [lenb, numb, numt]:
			### if(is.list(out$vbsc)){
			### 	out$vbsc = lapply(out$vbsc, function(x) if(length(dim(x)) == 2) array(x, dim=c(dim(x),1)) else x)
			### }
			### #else{
			### #	dim(out$vbsc) = c(length(out$vbsc)/500, 500)
			### #}
			### 		
			### # Check whether the dimensions of each time step are identical:
			### if(is.list(out$vbsc)){
			### 	out$vbsc = mergeListKeepDimensions(out$vbsc, pad=pad, split=split, add1=length(dim(out$vbsc[[1]])) == 2)
			### }
			### else if(length(out$vbsc) == 0){
			### 	warning(paste("Volume backscattering coefficient 'vbsc' not present for time step (s) ", paste(t, collapse=", "), sep=""))
			### }
			
			# Apply TVG:
			if(!TVG || TVG.exp != 2){
				if(!TVG){
					TVG.exp = 0
				}
				# Read beams-files and ctd-files:
				if(all(relevantbeamsvar %in% names(out))){
					beams = out[relevantbeamsvar]
				}
				else{
					beams = read.event_read_files(files=filelist, filesind=fInd$beams, tlist=tlist, var=beamsvar)
				}
				if(all(relevantctdvar %in% names(out))){
					ctd = out[relevantctdvar]
				}
				else{
					ctd = read.TSDs(filelist[fInd$ctd], var=relevantctdvar, msg=FALSE)
				}
				
				out$vbsc = apply.TVG(out$vbsc, beams=c(beams, ctd), rm=TRUE, TVG.exp=2)
				out$vbsc = apply.TVG(out$vbsc, beams=c(beams, ctd), rm=FALSE, TVG.exp=TVG.exp)
			}
				
			# Remove noise if required:
			if(any(!isTRUE(bgns), !isTRUE(pdns), !isTRUE(nrns), !isTRUE(hins))){
				# Read beams-files and ctd-files:
				if(all(relevantbeamsvar %in% names(out))){
					beams = out[relevantbeamsvar]
				}
				else{
					beams = read.event_read_files(files=filelist, filesind=fInd$beams, tlist=tlist, var=beamsvar)
				}
				if(all(relevantctdvar %in% names(out))){
					ctd = out[relevantctdvar]
				}
				else{
					ctd = read.TSDs(filelist[fInd$ctd], var=relevantctdvar, msg=FALSE)
				}
				# Get the noise variables:
				noisefiles = noise.path(event=event, esnm=esnm, dir.data=dir.data, utim=TIME$U000[t], noisevar=c("bgns","nrnp", "hins"))
				noise = NULL
				if(!isTRUE(bgns)){
					bgns = getbgns(if(is.character(bgns)) bgns else noisefiles)
					noise = c(noise, "bgns")
				}
				if(!isTRUE(pdns)){
					pdns = getpdns(if(is.character(pdns)) pdns else noisefiles, relevantpdnsvar=relevantpdnsvar)
					noise = c(noise, "pdns")
				}
				if(!isTRUE(nrns)){
					nrns = getnrns(if(is.character(nrns)) nrns else noisefiles, mode=c("p", "a"))
					noise = c(noise, "nrns")
				}
				if(!isTRUE(hins)){
					hins = gethins(if(is.character(hins)) hins else noisefiles)
					noise = c(noise, "hins")
				}
				# Calculate the noise:
				out$tlns = read.event_generate_noise(c(list(vbsc=out$vbsc), beams, ctd, bgns, pdns, hins, nrns), noise=noise, t=t, nsind=nsind, cruise=cruise, esnm=esnm, dir.data=dir.data, hins_add=hins_add, phase=TRUE, pdns_scale=pdns_scale, TVG=TVG)
				# Subtract the noise:
				out$vbsc = out$vbsc-c(out$tlns)
			}
				
			# Smooth along beams, if required:
			if(length(kern)>0 && kern>0){
				if(is.integer(kern)){
					out$vbsc = medSmooth1(out$vbsc, kern, ...)
				}
				else{
					out$vbsc = kernSmooth1(out$vbsc, kern=kern)
				}
			}
				
			# Generate uniformly distributed points, if required:
			if(readpsr){
				# Define the list of variables used as input to pplot3d.TSD(), given in the order used in that function:
				thisl = list(data=out, var=c(var, oldvar))
				# Define the variables present in '...' but not in the list 'thisl':
				otherl = ll[setdiff(names(ll), names(thisl))]
				# Set nlim to Inf:
				otherl$nlim = Inf
				thispsr = do.call("pplot3d.TSD", c(thisl, otherl))
				# Add to the output:
				sumatlevels = sapply(lapply(dim_all(thispsr$psxr), unlist), sum)
				if(all(sumatlevels == 0)){
					warning("The plotting region contains no regenerated data (possibly 'range' discards all the data, or 'acca' or 'N' is set too high)")
				}
				out$psxr = thispsr$psxr[[1]]
				thispsr$psxr = NULL
				out$psyr = thispsr$psyr[[1]]
				thispsr$psyr = NULL
				out$pszr = thispsr$pszr[[1]]
				thispsr$pszr = NULL
			}
				
			# Generate sv values in a cartesian grid:
			if(any(c("vbsC", "vbsA") %in% var)){
				#numt = length(out$psxr[[1]])
				####################### HERE THE out$psxr is a simple vector, and numt = length(out$psxr) is WRONG #################
				numt = length(out$psxr)
				out[c("psxg", "psyg", "pszg", "nbrp", "vbsC")] = rep(list(vector("list", length(numt))), 5)
				#out$psxg = vector("list", length(numt))
				#out$psyg = vector("list", length(numt))
				#out$pszg = vector("list", length(numt))
				#out$nbrp = vector("list", length(numt))
				#out$vbsC = vector("list", length(numt))
				out$ggsz = ggsz
				
				# Run through the time steps:
				for(i in seq_len(numt)){
					# Get the range of psxr, psyr, and pszr:
					# Midpoints:
					midGgridx = seq(floor(min(out$psxr[[i]])/ggsz), ceiling(max(out$psxr[[i]])/ggsz))*ggsz
					midGgridy = seq(floor(min(out$psyr[[i]])/ggsz), ceiling(max(out$psyr[[i]])/ggsz))*ggsz
					midGgridz = seq(floor(min(out$pszr[[i]])/ggsz), ceiling(max(out$pszr[[i]])/ggsz))*ggsz
					# Expand the midpoints:
					midGdrid = as.matrix(expand.grid(midGgridx, midGgridy, midGgridz))
					out$psxg[[i]] = midGdrid[, 1]
					out$psyg[[i]] = midGdrid[, 2]
					out$pszg[[i]] = midGdrid[, 3]
					# Create the grid:
					Ggridx = c(midGgridx[1]-ggsz/2, midGgridx + ggsz/2)
					Ggridy = c(midGgridy[1]-ggsz/2, midGgridy + ggsz/2)
					Ggridz = c(midGgridz[1]-ggsz/2, midGgridz + ggsz/2)
					# Find intervals:
					inGgridx = findInterval(out$psxr[[i]], Ggridx)
					inGgridy = findInterval(out$psyr[[i]], Ggridy)
					inGgridz = findInterval(out$pszr[[i]], Ggridz)
					# Count the points:
					flatgridindices = arr.ind2ind(cbind(inGgridx, inGgridy, inGgridz), c(length(midGgridx), length(midGgridy), length(midGgridz)))
					counts = tabulate(flatgridindices,length(out$psxg[[i]]))
					counts[counts == 0] = NA
					# Insert into the count variable:
					out$nbrp[[i]] = counts
					# Convert back to sv:
					out$vbsC[[i]] = out$nbrp[[i]]*thispsr$finalacca[i]/(ggsz^3)
				}
				
				if("vbsA" %in% var){
					if(all(c("psxx", "psyx", "pszx") %in% names(ll$range))){
						dimnamesx = seq(min(ll$range$psxx), max(ll$range$psxx), ggsz)
						dimnamesy = seq(min(ll$range$psyx), max(ll$range$psyx), ggsz)
						dimnamesz = seq(min(ll$range$pszx), max(ll$range$pszx), ggsz)
					}
					else{
						dimnamesx = seq(floor(min(unlist(lapply(out$psxg,min))) / ggsz) * ggsz, ceiling(max(unlist(lapply(out$psxg,max))) / ggsz) * ggsz, ggsz)
						dimnamesy = seq(floor(min(unlist(lapply(out$psyg,min))) / ggsz) * ggsz, ceiling(max(unlist(lapply(out$psyg,max))) / ggsz) * ggsz, ggsz)
						dimnamesz = seq(floor(min(unlist(lapply(out$pszg,min))) / ggsz) * ggsz, ceiling(max(unlist(lapply(out$pszg,max))) / ggsz) * ggsz, ggsz)
					}
					out$vbsA = zeros(length(dimnamesx), length(dimnamesy), length(dimnamesz), numt)
					dimnames(out$vbsA) = list(dimnamesx, dimnamesy, dimnamesz, t)
						
					# Run through the time steps:
					for(i in seq_len(numt)){
						# Discard points outside of the range:
						validgridcells = out$psxg[[i]] >= min(dimnamesx) & out$psxg[[i]] <= max(dimnamesx)
						validgridcells = validgridcells & out$psyg[[i]] >= min(dimnamesy) & out$psyg[[i]] <= max(dimnamesy)
						validgridcells = validgridcells & out$pszg[[i]] >= min(dimnamesz) & out$pszg[[i]] <= max(dimnamesz)
						indx = (out$psxg[[i]][validgridcells]-min(dimnamesx))/ggsz + 1
						indy = (out$psyg[[i]][validgridcells]-min(dimnamesy))/ggsz + 1
						indz = (out$pszg[[i]][validgridcells]-min(dimnamesz))/ggsz + 1
						out$vbsA[cbind(indx, indy, indz, i)] = out$vbsC[[i]][validgridcells]
					}
				}
			}
			# Calibration factor:
			if(cal != 1){
				if(is.list(out$vbsc)){
					out$vbsc = lapply(out$vbsc, function(x) x*c(matrix(cal,nrow=nrow(x), ncol=ncol(x))))
				}
				else{
					out$vbsc = out$vbsc*c(matrix(cal,nrow=nrow(out$vbsc), ncol=ncol(out$vbsc)))
				}
			}
				
			# Prepare a data frame of only the segmented data:
			if(isTRUE(mask) && (length(out$vxIX) || length(out$sgsc))){
				# Declare the mask:
				mask <- list()
				
				# Extract the segmentation indices:
				ind <- if(length(out$sgsc)) out$sgsc else out$vxIX
				if(!is.list(ind)){
					ind <- list(ind)
				}
				# Get the number of values per time step:
				numtPerPing <- sapply(ind, length)
				# Get the time:
				if(length(out$utim)){
					datetime <- as.POSIXct(out$utim, origin="2970-01-01")
					date <- format(datetime, "%Y-%m-%d")
					time <- format(datetime, "%H-%H-%OS3")
					mask$date <- rep(date, numtPerPing)
					mask$time <- rep(time, numtPerPing)
					mask$utim <- rep(out$utim, numtPerPing)
				}
				if(length(out$indt)){
					mask$indt <- rep(out$indt, numtPerPing)
				}
				
				# Get array indices:
				if(length(out$lenb) || length(out$numb)){
					l <- lapply(ind, ind2arr.ind, c(max(out$lenb), max(out$numb)))
					arr.ind <- do.call("rbind", l)
					mask$indb <- arr.ind[,2]
					mask$indr <- arr.ind[,1]
				}
				
				# Position data:
				if(length(out$psxx)){
					if(length(dim(out$psxx))==2){
						dim3 <- c(dim(out$psxx), 1)
						dim(out$psxx) <- dim3
						dim(out$psyx) <- dim3
						dim(out$pszx) <- dim3
					}
					mask$psxx <- unlist(lapply(seq_along(ind), function(i) out$psxx[,,i][ind[[i]]]), use.names=FALSE)
					mask$psyx <- unlist(lapply(seq_along(ind), function(i) out$psyx[,,i][ind[[i]]]), use.names=FALSE)
					mask$pszx <- unlist(lapply(seq_along(ind), function(i) out$pszx[,,i][ind[[i]]]), use.names=FALSE)
				}
				
				# Acoustic data:
				if(length(out$vbsc)){
					if(length(dim(out$vbsc))==2){
						dim(out$vbsc) <- c(dim(out$vbsc), 1)
					}
					mask$vbsc <- unlist(lapply(seq_along(ind), function(i) out$vbsc[,,i][ind[[i]]]), use.names=FALSE)
				}
				
			out$mask <- as.data.frame(mask)
			}
			
			# Extracting a subset if required:
			if(length(ll$ind)>0){
				ll$ind = ind.expand(ll$ind, dim_all(out$vbsc))
				out$vbsc = out$vbsc[ll$ind[[1]], ll$ind[[2]], , drop=FALSE]
			}
				
			# Adding Sv-values to the output list, if required:
			if("mvbs" %in% var){
				out$mvbs = 10*log10(out$vbsc)
			}
			# Removing the sv-values if required:
			if(!any(c("vbsc", "pings") %in% var)){
				out$vbsc = NULL
			}
		} # End of readsv
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("vbsc", "mvbs", "pings", "psxr", "psyr", "pszr", "psxg", "psyg", "pszg", "nbrp", "vbsC", "vbsA"))
	}	
	##################################
	##################################
	
	
	####################################
	### (8) Reading electric angles: ###
	####################################
	# Read the relevant pings files to obtain acoustic data. The method in the following block may be changed to a simpler method in the future, inspired by how the school-files are read below:
	if(any(c("angl", "angt") %in% var)){
		# Only read if 't' is not empty:
		if(length(t)>0){
			suppressWarnings(out$angl <- read.TSDs(filelist[fInd$pings], t=tlist[fInd$pings], var="angl", dimension=TRUE, merge=merge, indt=FALSE, msg=FALSE)$angl)			
			suppressWarnings(out$angt <- read.TSDs(filelist[fInd$pings], t=tlist[fInd$pings], var="angt", dimension=TRUE, merge=merge, indt=FALSE, msg=FALSE)$angt)
			
			# Check whether the dimensions of each time step are identical:
			if(is.list(out$angl)){
				out$angl = mergeListKeepDimensions(out$angl, pad=pad, split=split, add1=length(dim(out$angl[[1]])) == 2)
			}
			else if(length(out$angl) == 0){
				warning(paste("Volume backscattering coefficient 'angl' not present for time step (s) ", paste(t, collapse=", "), sep=""))
			}
			# Check whether the dimensions of each time step are identical:
			if(is.list(out$angt)){
				out$angt = mergeListKeepDimensions(out$angt, pad=pad, split=split, add1=length(dim(out$angt[[1]])) == 2)
			}
			else if(length(out$angt) == 0){
				warning(paste("Volume backscattering coefficient 'angt' not present for time step (s) ", paste(t, collapse=", "), sep=""))
			}
		}
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c("angl", "angt"))
	}
	####################################
	####################################

	
	#####################################
	### (9) Reading school variables: ###
	#####################################
	# Reading school files:
	if(any(var %in% c(staticschoolnames, dynschoolnames))){
		# Obtaining school file types:
		schooltype = echoIBM.getSchoolfileType(filelist[fInd$school], dynschoolnames, staticschoolnames)
		fInd$dynschool = schooltype$schooltypeD==1
		fInd$staticschool = schooltype$schooltypeS==1
		
		# Read static school variables:
		suppressWarnings(out <- c(out, read.TSDs(filelist[fInd$staticschool], var=var, t=1, msg=FALSE)))
			
		# Read dynamic school variables:
		if(any(dynschoolnames %in% var) && length(filelist[fInd$dynschool])>0){
			# Adding school-values to the output list:
			thesevar = intersect(var, dynschoolnames)
			out[thesevar] = rep(list(vector("list", length(t))), length(thesevar))
			# Read the school-files using the time step list 'tlist', setting dimension to TRUE to ensure that the data are preserved in time steps, clean = FALSE to ensure that all data is read (as clean = TRUE removes variables read from more than the first file if the same variable is present in more than one file), merge = TRUE to merge the time steps togeather for each variable, and indt = FALSE because of the way 'tlist' is defined:
			suppressWarnings(thesepings <- read.TSDs(filelist[fInd$dynschool], t=tlist[fInd$dynschool], var=var, dimension=TRUE, merge=merge, indt=FALSE, msg=FALSE))
			# Remove variables not present in the data and add the ones that are present to the output:
			out[names(thesepings)] = thesepings
			out[setdiff(thesevar, names(thesepings))] = NULL
			
			# If 'rtxf' or 'rtzf' are present in 'var' but not in 'out', an attempt is made to extract these rotation angles from the velocity vectors:
			if(!any(is.null(out$vlxf), is.null(out$vlyf), is.null(out$vlzf))){
				out = vl2rt.TSD(out, var=var)
			}
			else if(("rtxf" %in% var && is.null(out$rtxf)) || ("rtzf" %in% var && is.null(out$rtzf))){
				suppressWarnings(outforrot <- read.TSDs(filelist[fInd$dynschool], t=tlist[fInd$dynschool], var=c("vlxf", "vlyf", "vlzf"), dimension=TRUE, clean=FALSE, merge=merge, indt=FALSE, msg=FALSE))
				if(any(is.null(outforrot$vlxf), is.null(outforrot$vlyf), is.null(outforrot$vlzf))){
					warning("Fish rotation angles not available")
				}
				else{
					temp = vl2rt.TSD(outforrot, var=var)
					#out[names(temp)] = temp
					addnames <- setdiff(names(temp), names(out))
					out[addnames] = temp[addnames]
				}
			}
		}
		# Remove from 'var' the variables that have been read:
		var = setdiff(var, c(staticschoolnames, dynschoolnames))
	}
	#####################################
	#####################################
	
		
	#####################################
	### (10) Reading other variables: ###
	#####################################
	# Read all the variables that are not in .pings-files, and that are not school variables:
	# merge = merge is needed to comply with the input, and indt = FALSE is a consequence of the way 'tlist' is defined (read.event() always considers general time step indexes and 'tlist' is defined according to this requirement):
	if(other && length(var)>0){
		# Read the other variables:
		otherout = read.TSDs(filelist[fInd$notpings], t=tlist[fInd$notpings], var=var, dimension=TRUE, merge=merge, clean=FALSE, indt=FALSE, msg=FALSE, addNvar=FALSE)
		
		# Add empty list elements where no data were present:
		if(merge){
			otherout = read.event_insertNA(otherout, TIME, fInd$notpings, filelist, t, tlist, var, var)
		}
		
		# Add to the output:
		out[names(otherout)] = otherout
	}
	#####################################
	#####################################
	
	# Transform from decibar to Pascal:
	if(Paout){
		if(!is.null(out$ihpr)){
			out$ihpr = out$ihpr*10000
		}
		if(!is.null(out$hpr0)){
			out$hpr0 = out$hpr0*10000
		}
	}
	
	# Drop redundant dimensions of the elements of 'out' if required:
	if(drop.out){
		#out = lapply(out, drop)
		out = lapply(out, function(x) if(is.list(x) && length(x) == 1) unlist(x) else drop(x))
	}
			
	# Stripping the duplicated elements:
	if(TIMEout){
		if(strip.out){
			out = out[!duplicated(names(out))]
			out["TIME"] = list(TIME)
		}
		else{
			out["TIME"] = list(TIME)
		}
	}
	else{
		if(strip.out){
			out = out[!duplicated(names(out))]
		}
		else{
			out = out
		}
	}
	
	# Add the variables given in 'adds':
	if(length(adds)>0){
		out[names(adds)] = adds
	}

	# Add info and return:
	if(length(out)){
		out$info = info.TSD(names(out))
	}
	out
	##################################################
	##################################################
}
