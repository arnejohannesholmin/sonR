#*********************************************
#*********************************************
#' Write data to a netCDF4 file
#'
#' \code{writeNetCDF4} Writes a list of data to a netCDF4 file.
#' \code{tsd2NetCDF4} Writes a list of data to a netCDF4 file.
#' \code{pings2NetCDF4} Convenience function for converting .pings files to netCDF4 files.
#' \code{seg2NetCDF4} Convenience function for converting .seg files to netCDF4 files.
#'
#' @param data			The list of data to write to the netCDF4 file.
#' @param con			The path to the netCDF4 file.
#' @param units			An optional list of units, named by the names of \code{data}.
#' @param dims			An optional list of dimensions, named by the names of \code{data}, where each element must be a list of the elements Name and Unit as vectors of the names and units of the dimensions of the element of \code{data} that the dimensions correspond to (example: dims = list(vbsc = list(Name = c("sample", "beam", "ping"), Unit = c("SampleIndex", "BeamIndex", "PingIndex")))).
#' @param tsdFile		The path to a TSD file.
#' @param netCDF4dir	The path to the directory in which to put a netCDF4 file.
#' @param event			The path to the event holding TSD files.
#'
#' @return
#'
#' @importFrom ncdf4 ncdim_def ncvar_def nc_create ncvar_def ncvar_put nc_close
#'
#' @rdname writeNetCDF4
#' @export
#'
writeNetCDF4 <- function(data, con, units = NULL, dims = NULL) {
	
	# Function to put one data element to netCDF4:
	putData <- function(ind, vars, data, nc) {
		ncdf4::ncvar_put(
			nc = ncnew, 
			varid = vars[[ind]], 
			vals = data[[ind]]
		)
	}
	# Function fo extract the given names and units of the dimensions of a specific element of the data:
	getDimsDFOne <- function(name, dims, data) {
		if(name %in% names(data)) {
			data.frame(
				Dim = dim(data[[name]]), 
				Name = dims[[name]]$Name, 
				Unit = dims[[name]]$Unit, 
				stringsAsFactors = FALSE
			)
		}
		else {
			NULL
		}
	}
	
	# Get the dimensions of the variable:
	if(length(dims) && is.list(dims)) {
		# Get the individual dimsDFs and rbind:
		dimsDF <- lapply(names(dims), getDimsDFOne, dims = dims, data = data)
		dimsDF <- do.call(rbind, dimsDF)
	}
	else {
		dimsDF <- NULL
	}
	dims <- defineNetCDF4Dims(data, dimsDF = dimsDF)
	
	# Define the variables (with dimensions):
	if(length(units) == 0) {
		units <- as.list(rep("NA", length(data)))
		#units <- as.list(rep(NA, length(data)))
		names(units) <- names(data)
	}
	vars <- lapply(
		names(data), 
		defineNetCDF4OneVar, 
		data = data, 
		dims = dims, 
		units = units
	)
	
	# Create the file:
	dir.create(dirname(con), showWarnings = FALSE, recursive = TRUE)
	ncnew <- ncdf4::nc_create(con, vars = vars)
	
	# Add the netCDF4 data:
	lapply(
		seq_along(data), 
		putData, 
		vars = vars, 
		data = data, 
		nc = ncnew
	)
	
	# Close the file:
	ncdf4::nc_close(ncnew)

	con
}
#'
#' @rdname writeNetCDF4
#' @export
#'
tsd2NetCDF4 <- function(tsdFile, netCDF4dir, units = NULL, dims = NULL) {
	
	# Read the vbsc data:
	data <- read.TSD(tsdFile, var = "all", t = "all")
	
	# Get the path to the netCDF4 file:
	netCDF4file <- paste0(basename(tools::file_path_sans_ext(tsdFile)), "_", tools::file_ext(tsdFile), ".nc")
	netCDF4filePath <- file.path(netCDF4dir, netCDF4file)
	
	# Write the file:
	writeNetCDF4(
		data = data, 
		con = netCDF4filePath, 
		units = units, 
		dims = dims
	)
}
#'
#' @rdname writeNetCDF4
#' @export
#'
pings2NetCDF4 <- function(event) {
	
	# Get the event:
	if(is.list(event)) {
		event <- event$path
	}
	
	# Get the pings files:
	pingsFiles <- list.files(event, full.names = TRUE, pattern = "\\.pings$")
	
	# Define the dimensions:
	dims <- dims_vbsc
	
	# Get and create the netCDF4 directory:
	netCDF4dir <- file.path(dirname(event), "netCDF4")
	dir.create(netCDF4dir, showWarnings = FALSE, recursive = TRUE)
	
	
	sapply(pingsFiles, tsd2NetCDF4, netCDF4dir = netCDF4dir, dims = dims)
}
#'
#' @rdname writeNetCDF4
#' @export
#'
seg2NetCDF4 <- function(event) {
	
	# Get the event:
	if(is.list(event)) {
		event <- event$path
	}
	
	# Get the pings files:
	segFiles <- list.files(event, full.names = TRUE, pattern = "\\.seg$")
	
	# Define the dimensions:
	dims <- dims_segM
	
	# Get and create the netCDF4 directory:
	netCDF4dir <- file.path(dirname(event), "netCDF4")
	dir.create(netCDF4dir, showWarnings = FALSE, recursive = TRUE)
	
	
	sapply(segFiles, tsd2NetCDF4, netCDF4dir = netCDF4dir, dims = dims)
}


dims_vbsc <- list(
	vbsc = list(
		Name = c("sample", "beam", "ping"), 
		Unit = c("SampleIndex", "BeamIndex", "PingIndex")
	)
)
dims_segM <- list(
	segM = list(
		Name = c("sample", "beam", "ping"), 
		Unit = c("SampleIndex", "BeamIndex", "PingIndex")
	)
)


defineNetCDF4Dims <- function(data, dimsDF = NULL) {
	
	# Get the dimensions of the variable:
	allDims <- unique(unlist(lapply(data, dim_all)))
	allDimSeq <- lapply(allDims, seq_len)
	
	# Get the dimension names and units from the input dimsDF data.frame (with columns Dim, Name, Unit):
	if(length(dimsDF)) {
		m <- match(dimsDF$Dim, allDims)
		dimsNames <- dimsDF$Name[m]
		dimsUnits <- dimsDF$Unit[m]
	}
	# ... or apply default values:
	else {
		dimsNames <- paste0("Dim", seq_along(allDims))
		dimsUnits <- paste0("Unit", seq_along(allDims))
	}
	
	# Define the dimensions of the data:
	dims <- mapply(ncdf4::ncdim_def, name = dimsNames, units = dimsUnits, vals = allDimSeq, SIMPLIFY = FALSE)
	
	dims
}

defineNetCDF4OneVar <- function(name, data, dims, units = NULL) {
	
	if(length(units) == 0) {
		units <- "Missing"
	}
	else if(is.list(units)){
		units <- units[[name]]
	}
	
	# Recognize the dimensions of each element of the data list:
	thisdim <- dim_all(data[[name]])
	allDims <- sapply(dims, "[[", "len")
	m <- match(thisdim, allDims)
	
	# Define the vbsc netCDF4 data:
	var <- ncdf4::ncvar_def(
		name = name, 
		units = units, 
		dim = dims[m]
	) 
	
	var
}
