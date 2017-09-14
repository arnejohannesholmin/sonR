#*********************************************
#*********************************************
#' (Internal) Get the name of the acoustic system form the raw data.
#'
#' @return
#'
#' @examples
#' \dontrun{}
#'
#' @export
#' @rdname EKRaw2TSD_getesnm
#' 
EKRaw2TSD_getesnm <- function(thisd, numt){
	if( length(grep("MBS", thisd$header$soundername) > 0) || length(grep("MS70", thisd$header$soundername) > 0) || thisd$header$transceivercount == 500 ){
		rep("MS70",numt)
		}
	else if( length(grep("MBES", thisd$header$soundername) > 0) || length(grep("ME70", thisd$header$soundername) > 0) ){
		rep("ME70",numt)
		}
	else if( length(grep("ER60", thisd$header$soundername) > 0) || length(grep("EK60", thisd$header$soundername) > 0) ){
		rep("EK60",numt)
		}
	else if( length(grep("SX80", thisd$header$soundername) > 0) ){
		rep("SX80",numt)
		}
	else if( length(grep("SH80", thisd$header$soundername) > 0) ){
		rep("SH80",numt)
		}
	else if( length(grep("SU80", thisd$header$soundername) > 0) ){
		rep("SU80",numt)
		}
	else if( length(grep("SX90", thisd$header$soundername) > 0) ){
		rep("SX90",numt)
		}
	else if( length(grep("SH90", thisd$header$soundername) > 0) ){
		rep("SH90",numt)
		}
	else if( length(grep("SU90", thisd$header$soundername) > 0) ){
		rep("SU90",numt)
		}
	else{
		rep(thisd$header$soundername,numt)
		}
	}
