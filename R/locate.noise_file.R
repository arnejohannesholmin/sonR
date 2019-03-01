#*********************************************
#*********************************************
#' (Internal) Locate noise file.
#'
#' @param event		The path to the event.
#' @param noisevar	A vector of noise variabes to find.
#' @param corvar	A vector of noise correlation variabes to find.
#' @param esnm		The name of the acoustic instrument for which to find the noise files.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @importFrom TSD read.TSDs strff
#' @importFrom utils tail
#'
#' @export
#' @rdname locate.noise_file
#' 
# Function for locating noise files in an event or directory:
locate.noise_file <- function(event, noisevar=c("bgns","nrnp"), corvar=c("crb1","olpb"), esnm="MS70"){
	if(length(event)>0 && file.exists(event)){
		# If a UNIX_time file exists, read this, otherwise generate one:
		TIME <- UNIX_time(event, file=TRUE)
		if(length(TIME) == 0){
			return(list())
			}
		readesnm <- read.TSDs(unlist(TIME$f000), var="esnm", clean=FALSE, merge=FALSE, info=FALSE)
		readesnm <- sapply(readesnm, function(x) strff(esnm[1], x) | length(x)==0)
		
		# Match against the relevant variables:
		arenoise <- unlist(lapply(TIME$l000, function(x) any(noisevar %in% x))) & readesnm
		arecor <- unlist(lapply(TIME$l000, function(x) any(corvar %in% x))) & readesnm
		
		if(sum(arenoise) || sum(arecor)){
			list(noisefiles=TIME$f000[arenoise], corfiles=TIME$f000[arecor], noise_utim=sapply(TIME$u000[arenoise],tail,1), cor_utim=sapply(TIME$u000[arecor],tail,1))
			}
		else{
			list()
			}
		}
	else{
		list()
		}
	}
