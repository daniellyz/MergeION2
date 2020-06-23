#' Generating consensus spectra
#'
#' Function that merges a list of normalized spectra that belong to the same ID
#'
#' @param splist List of extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity
#' @param mz_window  m/z tolerance window (in da) for spectra alignment
#'
#' @return
#' \itemize{
#'   \item{new_spectrum:}{The aligned and averaged consensus spectra: data matrix with two columns (averaged m/z and intensity)}
#'   \item{I_matrix:}{The aligned intensity matrix: intensity of aligned m/z in in each MS record}
#'   \item{dev_matrix:}{The mass deviation matrix: mass deviation of each m/z in each MS record}
#' }
#' 
#' @export
#'
average_spectrum<-function(splist, mz_window = 0.01){
  
  if (!is.list(splist)){
    stop("Please provide a list of mass spectra!")
  }
  
  if (length(splist)<2){
    stop("At least 2 spectra is needed for consensus spectra generation!")
  }
  
  ####################################################
  ### Combine all spectra into a m/z-ordered matrix:##
  ####################################################
  
  NBS = length(splist) # Number of spectra to average
  samples = c()
  for (i in 1:NBS){samples = rbind(samples,cbind(i,splist[[i]]))} # Create a big matrix
  samples = data.frame(samples)
  colnames(samples)=c("File_num","mz","I")
  samples = samples[order(samples$mz,decreasing = F),]
  
  ###############################################
  ### Find cluster of masses and assign labels:##
  ###############################################
  
  mz_labels = cut_mz_list(samples$mz, mz_window) # Cut according to mass
  mz_features = unique(mz_labels)
  N_feature=length(mz_features) # Number of molecular features
  
  #################################################
  ### Combine and average m/z with the same label##
  #################################################
  
  I_matrix = matrix(0,N_feature,NBS) # Intensity matrix N_feature x Nb of spectra
  avg_mzlist=rep(0,N_feature) # Averaged mass list
  dev_matrix = matrix(0,N_feature,NBS)
  
  for (i in 1:N_feature){
    
    valid = which(mz_labels == mz_features[i])
    avg_mzlist[i] = mean(samples$mz[valid])
    dev_matrix[i, samples$File_num[valid]] = samples$mz[valid] - avg_mzlist[i]
    
    for (j in 1:length(valid)){ # Actualize intensity matrix
      I_matrix[i,samples$File_num[valid[j]]]= I_matrix[i,samples$File_num[valid[j]]]+samples$I[valid[j]]
  }}
  
  ########################
  ### Output new spectrum:
  ########################
  
  new_intensity = apply(I_matrix, 1, mean)
  new_spectrum = cbind(avg_mzlist,new_intensity)
  colnames(new_spectrum)=NULL
  rownames(new_spectrum)=NULL
  new_spectrum = new_spectrum[order(new_spectrum[,1]),]
  
  return(list(new_spectrum = new_spectrum, I_matrix = I_matrix, dev_matrix = dev_matrix))
}

##########################
### Internal functions:###
##########################

cut_mz_list<-function(mzlist, mz_window){
  
  N=length(mzlist)
  
  f=1
  mz_feature=c(0, N) 
  t0 = 1 # Start index of a cluster
  
  for (k in 2:N){
    min_mz = min(mzlist[t0:(k-1)])
    avg_mz = mean(mzlist[t0:(k-1)])
    
    if (mzlist[k] - min_mz > mz_window & mzlist[k] - avg_mz > mz_window/2){
      mz_feature[t0:(k-1)] = f 
      f = f + 1
      t0 = k
    }
  }
  mz_feature[t0:N] = f
  return(mz_feature)
}

custom.dist <- function(my.list, my.function) {
  n <- length(my.list)
  mat <- matrix(0, ncol = n, nrow = n)
  colnames(mat) <- rownames(mat) <- names(my.list)
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      mat[i,j] <- my.function(my.list[i],my.list[j])
    }}
  return(as.dist(mat))
}

