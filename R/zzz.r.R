# Project: MergeION2
# 
# Author: yzhan482
###############################################################################

#' @importFrom utils packageDescription
.onLoad <- function(libname, pkgname){
	options(error = NULL)
	packageStartupMessage(paste("\nMergeion version ", packageDescription("Mergeion")$Version, 
					"\n", sep = ""))
	
	# install github deps
	installedPkgs <- .packages(all.available = TRUE)
	if(! "dqshiny" %in% installedPkgs) remotes::install_github("daqana/dqshiny")
	if(! "compMS2Miner" %in% installedPkgs) remotes::install_github("WMBEdmands/compMS2Miner")
}

