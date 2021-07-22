#' Extracting or generating representative or consensus scans for each compound ID
#'
#' Function used by library_generator to generate consensus scans
#' 
#' @export
#'
#' 
process_consensus<-function(input_library, method = c("most_recent", "consensus", "consensus2", "common_peaks")[1], 
                            consensus_window = 0.01, relative=0.01, max_peaks = 200){

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  if (!(method  %in% c("consensus","consensus2","common_peaks","most_recent"))){
    stop("The library processing method does not exist!")
  }
  
  message("Generating consensus library...")
  
  ####################################
  ### Read and check input library:###
  ####################################
  
  input_library = library_reader(input_library)
  complete_library = input_library$complete
  
  spectrum_list = complete_library$sp
  metadata = complete_library$metadata

  ###############
  ### Grouping:##
  ###############
  
  labels = paste(metadata$IONMODE, metadata$MSLEVEL, sep="-")
  groups = unique(labels)
  
  #################
  ### Initialize:##
  #################
  
  new_spectrum_list = list()
  new_metadata = c()
  NN = 0
  
  #############################
  ### Combine scan per group:##
  #############################

  for (gg in groups){

    index1 = which(labels==gg)
    
    metadata1= metadata[index1,,drop=FALSE]
    spectrum_list1 = spectrum_list[index1]
    ID_list1=unique(metadata1$ID)
    
    for (ID in ID_list1){
      
      print(which(ID_list1==ID))
      selected_rows = which(metadata1$ID==ID)
      NSR = length(selected_rows)
      
      sub_metadata = metadata1[selected_rows,,drop=FALSE]
      sub_spectrum_list = spectrum_list1[selected_rows]
      
      # Representative scans in the library:

      wm = which.max(as.numeric(sub_metadata$SCANS)) 
      new_metadata = rbind.data.frame(new_metadata,sub_metadata[wm,,drop=FALSE]) # Update metadata
      NN = NN+1

      # Append spectra list if no need for spectra merge
      
      if (method == "most_recent" || NSR==1){
        new_spectrum_list[[NN]]=sub_spectrum_list[[wm]]
      }
      
      # Append consensus spectrum:
      
      if (method %in% c("consensus", "consensus2", "common_peaks") && NSR>1){
          output_consensus = average_spectrum(sub_spectrum_list, consensus_window)
          temp_spectrum = output_consensus$new_spectrum
          
          if (method == "common_peaks"){
            temp_zeros =  apply(output_consensus$I_matrix, 1, function(x) sum(x==0))
            temp_spectrum = temp_spectrum[temp_zeros==0,,drop=FALSE]
          }
          if (method == "consensus2"){
            temp_nz =  apply(output_consensus$I_matrix, 1, function(x) sum(x>0))
            temp_spectrum = temp_spectrum[temp_nz>=2,,drop=FALSE]
          }
          new_spectrum_list[[NN]] = temp_spectrum
      }
    }
  }

  #############################
  ### Denoising and Filtering:#
  #############################
  
  new_spectrum_list1 = list()
  included=c()
  n0=0
  
  for (i in 1:NN){
    
    sp0 = new_spectrum_list[[i]]
    sp1 = denoise_spectrum(sp0, max_peaks, relative)
    
    if (nrow(sp1)>1){
      included = c(included, i)
      n0 = n0 + 1
      new_spectrum_list1[[n0]]=sp1
    }
  }
  
  new_metadata = new_metadata[included,,drop=FALSE]
  new_metadata$PARAM_CONSENSUS = method
  
  ####################
  ### Return results:
  ####################

  consensus_library = list(metadata = new_metadata, sp = new_spectrum_list1)

  output_library = list(complete = complete_library, consensus = consensus_library, network = NULL)

  return(output_library)
}

#######################
### Internal function:#
#######################

# Keep top peaks and noramalzie

denoise_spectrum<-function(sp, max_peak, min_relative){
  
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

