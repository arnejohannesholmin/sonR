#*********************************************
#*********************************************
#' Merges time steps of acoustic data located in the neighbor directory "temp" of the "tsd"-directory given by 'event'.
#'
#' @param x			A list of data as returned from read.event, including vbsc.
#' @param filename	The path to the netcdf4 file.
#'
#' @return
#'
#' @examples
#' \dontrun{x <- read.event("~/Data/echoIBM/MachineLearning/Events/E021/EK60/tsd", var="vbsc", t = "all")}
#'
#' @importFrom ncdf4 ncdim_def ncvar_def nc_create ncvar_def ncvar_put nc_close
#'
#' @export
#' @rdname echoIBM
#'
vbsc2netcdf4 <- function(event) {
	
	if(is.list(event)) {
		event <- event$path
	}
	
	# Get the pings files:
	pingsFiles <- list.files(event, full.names = TRUE, pattern = "\\.pings$")
	
	netcdf4dir <- file.path(dirname(event), "netcdf4")
	dir.create(netcdf4dir, showWarnings = FALSE)
	
	sapply(pingsFiles, vbsc2netcdf4One, netcdf4dir = netcdf4dir)
}



vbsc2netcdf4One <- function(tsdFile, netcdf4dir) {
	
	# Read the vbsc data:
	vbsc <- read.TSD(tsdFile, var = "vbsc", t = "all")$vbsc
	
	netcdf4file <- paste(basename(tools::file_path_sans_ext(tsdFile)), "nc", sep = ".")
	netcdf4filePath <- file.path(netcdf4dir, netcdf4file)
	
	# Get the dimensions of the vbsc data:
	d <- dim(vbsc)
	dseq <- lapply(d, seq_len)
	names(dseq) <- c("sample", "beam", "ping")
	# Define the dimensions of the vbsc data:
	dims <- mapply(ncdf4::ncdim_def, name = names(dseq), units = "index", vals = dseq, SIMPLIFY = FALSE)
	
	# Define the vbsc netcdf4 data:
	var_temp <- ncdf4::ncvar_def(
		name = "vbsc", 
		units = "m^-1", 
		dim = dims
	) 
	
	# Create the file:
	ncnew <- ncdf4::nc_create(netcdf4filePath, list(var_temp))
	
	# Add the vbsc netcdf4 data:
	ncdf4::ncvar_put(
		nc = ncnew, 
		varid = var_temp, 
		vals = vbsc
	)
	
	# Close the file:
	ncdf4::nc_close(ncnew)

	netcdf4filePath
}
