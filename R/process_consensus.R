#' Extracting or generating representative or consensus scans for each compound ID
#' @inheritParams  library_generator
#' @param    IDsUpadated new compound ID which were updated in the \bold{complete} library yet to update consensus library for. 
#' 
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr group_by ungroup filter %>% mutate
process_consensus <- function(input_library, method = c("most_recent", "consensus", "consensus2", "common_peaks"), 
                            consensus_window = 0.01, relative = 0.01, max_peaks = 200, 
                            IDsUpdated = NULL){
  
  options(stringsAsFactors = FALSE)
  options(warn=-1)
  
  method <- match.arg(method)
  
  if (!(method  %in% c("consensus","consensus2","common_peaks","most_recent"))){
    stop("The library processing method does not exist!")
  }
  
  if(all(c("metadata", "sp") %in% names(input_library)))  input_library <- list(complete = input_library)
    
  if(is.null(input_library$complete)) stop(" Missing input library")
  
  missingInputConsensus <- is.null(input_library$consensus)
  
  # scenario 1: missing input consensus, and append_lib = NULL
  if(missingInputConsensus & is.null( IDsUpdated )){
    
    generateConcensusCondition <- TRUE
    message("Generating consensus library")
    
  }else if(missingInputConsensus & !is.null( IDsUpdated )){
    # scenario 2: missing input consensus, and append_lib != NULL
    #combine append_library to input_library
    
    generateConcensusCondition <- TRUE
    message("Missing consensus libary in input_library, generating consensus library for all compound present")
    
  }else if(!missingInputConsensus & is.null( IDsUpdated )){
    # scenario 3: input consensus is not missing, and append_lib = NULL
    generateConcensusCondition <- TRUE
    message("No additional compound IDs are provided, the consensus library will be updated using input_library")
  }else{
    # scenario 4: input consensus is not missing, and append_lib != NULL
    # the current consensus library will be udated
    generateConcensusCondition <- FALSE
    message("Updating existing concensus library") 
  }
  
  #############################################
  # make input_library a tibble of nested list
  ##############################################
  libraryTB <- tibble(sp = input_library$complete$sp, input_library$complete$metadata)
  #test <- libraryTB %>% filter(ID == "D024069") 
  #updateSpecOneID(listSpec =  test$sp, method = method, scans =  test$SCANS, consensus_window = 0.01)
  
  
  if( !generateConcensusCondition ) libraryTB <-   libraryTB %>% filter(ID %in%  IDsUpdated)
  
  concensusLibTemp <- libraryTB %>% group_by(ID, IONMODE,MSLEVEL) %>%
    mutate(spNew = list(updateSpecOneID(listSpec = sp, method = method, scans =  SCANS, consensus_window = consensus_window))) %>%
    filter(SCANS == max(SCANS)) %>%
    ungroup() 
  
  ##  denoise
  concensusLib <- list()
  denoisedSpectra <- sapply(concensusLibTemp$spNew, 
                            denoise_spectrum, max_peak = max_peaks, min_relative = relative, simplify = FALSE)
  validOnes <- !unlist(lapply( denoisedSpectra, is.null))
  concensusLib$sp <-   denoisedSpectra[validOnes] 
  concensusLib$metadata <-  concensusLibTemp[  validOnes, -1]
  # keep the original order
  concensusLib$metadata <-  concensusLib$metadata[, colnames(input_library$complete$metadata)]
  
  concensusLib$metadata$PARAM_CONSENSUS <- method
  ####################
  ### Return results:
  ####################
  
  if(generateConcensusCondition){
    
    output_library = list(complete = input_library$complete, consensus =   concensusLib, network = NULL)
    
    
  }else{
    
    inputConsensus <- input_library$consensus
    idx <-  inputConsensus$metadata$ID %in% IDsUpdated
    
    consensus_library = list(
                             sp = c(inputConsensus$sp[!idx],   concensusLib$sp),
                             metadata = rbind( inputConsensus$metadata[!idx,], concensusLib$metadata))
    
    output_library = list(complete = input_library$complete, consensus =  consensus_library, network = NULL)
    
    
  }
  
  
  return(output_library)
}

#######################
### Internal function:#
#######################

#' Keep top peaks and noramalzie
#' @inheritParams  library_generator
#' @param sp spectrum

denoise_spectrum <- function(sp, max_peak, min_relative){
  
  denoised_spectrum = matrix(c(0,0),1,2)
  
  if (nrow(sp)>0){
    
    # Filter top peaks:
    
    sp = sp[order(sp[,2], decreasing = T),,drop=FALSE]
    tops = min(max_peak, nrow(sp))  
    sp = sp[1:tops,,drop=FALSE]
    
    # Normalize to 100:
    
    sp1 = sp
    sp1[,2] = sp1[,2]/max(sp1[,2])*100
    
    # Relative Intensity filter:
    
    filter = which(sp1[,2]>=min_relative)
    sp = sp1  
    sp = sp[filter,,drop=FALSE]
    
    # Check validity:
    
    if (nrow(sp)>=2){
      sp = sp[order(sp[,1]),]
      sp[,2] = sp[,2]/max(sp[,2])
      denoised_spectrum = sp
    }
  }
  return(denoised_spectrum)
}


#' Genereate consensus library for a given compound ID
#' @param listSpec list of spectra
#' @inheritParams  library_generator
#' @param scans the "SCANS" from the metadata
#' ' @importFrom tibble as_tibble


updateSpecOneID <- function(listSpec, method = c("most_recent", "consensus", "consensus2", "common_peaks"),
                            scans, consensus_window){
  
  method <- match.arg(method)
  nrRecord <- length(scans)
  maxRecord <- which.max(scans)
  
  
  if (method == "most_recent" || nrRecord == 1){
    newSpectrum <- listSpec[[  maxRecord ]]
	newSpectrum <- newSpectrum[order(newSpectrum[,1]),,drop = FALSE]
  }else{
    
    if (method %in% c("consensus", "consensus2", "common_peaks") && nrRecord > 1){
      output_consensus <- average_spectrum(listSpec, consensus_window)
      newSpectrum <- output_consensus$new_spectrum
      
      if (method == "common_peaks"){
        temp_zeros <-  apply(output_consensus$I_matrix, 1, function(x) sum(x==0))
        newSpectrum <-   newSpectrum[temp_zeros==0,,drop=FALSE]
      }else if (method == "consensus2"){
        temp_nz <-  apply(output_consensus$I_matrix, 1, function(x) sum(x>0))
        newSpectrum <- newSpectrum[temp_nz>=2,,drop=FALSE]
      }
    }
  }
  newSpectrum <-  as_tibble(newSpectrum)
  colnames(newSpectrum)  <- c("m/z","Intensity")
  
  return(newSpectrum)
  
}

