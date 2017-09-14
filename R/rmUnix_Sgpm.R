#*********************************************
#*********************************************
#' Remove the files with basenames UNIX_time.tsd and sgPM.tsd. Used in UNIX_time() and read.event().
#'
#' @param x  A vector of file names.
#' @export
#' @rdname rmUNIX_SgPM
#'
rmUnix_Sgpm <- function(x){
	x[!tolower(basename(x)) %in% c("unix_time.tsd","sgpm.tsd")]
}