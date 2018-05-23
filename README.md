sonR R package
=====

This package contains many functions used by the packages echoIBM and cpplot3d.

Version: 1.2
Required R version: 3.3.3

Installation
=====

``` r
# Install the packages that sonR depends on. Note that this updates all the specified packages to the latest (binary) version. To skip installing already installed packages, run install.packages(setdiff(dep.pck, installed.packages()[,"Package"]), repos="http://cran.us.r-project.org") instead:
dep.pck <- c("devtools", "ccaPP", "data.table", "fpc", "pbapply", "SoDA", "XML")
install.packages(dep.pck, repos="http://cran.us.r-project.org")

# Install sonR and also the packages that sonR depends on which are on GitHub (by Holmin):
# On Windows you will need Rtools to complete the installations. Check if you have this by running Sys.getenv('PATH'), and go to https://cran.r-project.org/bin/windows/Rtools/ to install Rtools if not. Note that if you need to run R as administrator due to security settings, it is advised to install the pakcages in plain R, and not using Rstudio. Close Rstudio, open R and run the installation, and reopen Rstudio.

dep.pck.git <- c("arnejohannesholmin/TSD", "arnejohannesholmin/SimradRaw", "arnejohannesholmin/sonR")
# If you want to install the lastest development versions, run devtools::install_github(dep.pck.git, ref="develop") instead:
devtools::install_github(dep.pck.git)

```

# For changes log see https://github.com/arnejohannesholmin/sonR/NEWS

Examples
=====

``` r
# This is simply a demonstration of how the sonR package converts raw data to TSD files. The raw files are put in an event, which is a directory path CruiseName/"Events"/EventName/EchosounderName/"raw such as S2017114/Events/Event001/EK60/raw. Using EKRaw2TSD() the raw files are converted to TSD files in a directory "tsd" alongside the "raw" directory. We do this in the tempdir() here, but the used should select a directory in which to put the cruises using the function sonR::Acoustics_datasets_directory(INSERT_PATH_TO_ACOUSTIC_DATA_HERE).

# Set the directory of the acoustic data, here simply as the tempdir() but preferably another location:
dir <- tempdir()

ev <- generate.event(event="Event1", cruise="Cruise1", esnm="EK60", dir.type = c("raw", "tsd"), dir.data=dir)
evRaw <- ev[1]
evTSD <- ev[2]
# Add one raw file to the event:
echoSounderFile <- file.path(system.file("extdata", package="sonR"), "RedSlip-D20160915-T120914.raw")
file.copy(echoSounderFile, evRaw)
# Generate the TSD files, which are faster to read with R:
EKRaw2TSD(evRaw)


# Read Sv (mvbs) and sv (vbsc) data:
d <- read.event(evTSD, t=1:6, var=c("mvbs", "vbsc", "voxels"))
# Get description of frequently used variables:
info.TSD(labl.TSD("a")) # Acoustic variables
info.TSD(labl.TSD("vx")) # Voxel variables (position, volume)
info.TSD(labl.TSD("v")) # Vessel variables
info.TSD(labl.TSD("rb")) # (Relevant) beams variables

# These are arrays with dimension [samples, beams, pings]:
str(d)
hist(d$mvbs)
# Plot one ping of the first frequency:
plot(d$mvbs[,1], type="l")
# The same along depth:
plot(d$mvbs[,1], d$pszx[,1], type="l")

summary(c(d$mvbs))
summary(c(d$vbsc))
```

License
=====

The sonR package is licensed under the LGPL-3.)

