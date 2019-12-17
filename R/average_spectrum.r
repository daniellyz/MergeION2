#' Generating consensus spectra
#'
#' Internal function that merges a list of normalized spectra that belong to the same ID
#'
#' @param splist List of extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity
#' @param ppm_window  m/z tolerance window (in ppm) for spectra alignment
#' @param clean Logical. TRUE if peaks not present in all spectra are filtered.
#'
#' @return
#' new_spectrum ~ The aligned, (filtered) and averaged consensus spectra: data matrix with two columns (m/z and intensity)
#'
#' @export
#'
#' @importFrom stats cutree hclust as.dist

average_spectrum<-function(splist, ppm_window = 10, clean = F){

  if (!is.list(splist)){
    stop("Please provide a list of mass spectra!")
  }

  if (length(splist)<2){
    stop("At least 2 spectra is needed for consensus spectra generation!")
  }

  ####################################################
  ### Combine all spectra into a m/z-ordered matrix:
  ####################################################

  NBS = length(splist) # Number of spectra to average
  samples = c()
  for (i in 1:NBS){samples = rbind(samples,cbind(i,splist[[i]]))} # Create a big matrix
  samples = data.frame(samples)
  colnames(samples)=c("File_num","mz","I")
  samples = samples[order(samples$mz,decreasing = F),]

  ##############################################
  ### Find cluster of masses and assign labels:
  ##############################################

  mz_labels = cut_mz_list(samples$mz,ppm_window) # Cut according to mass
  mz_features = unique(mz_labels)
  N_feature=length(mz_features) # Number of molecular features

  ###############################################
  ### Combine and average m/z with the same label:
  ###############################################

  I_matrix=matrix(0,N_feature,NBS) # Intensity matrix N_feature x Nb of spectra
  avg_mzlist=rep(0,N_feature) # Averaged mass list

  for (i in 1:N_feature){

    valid = which(mz_labels == mz_features[i])
    avg_mzlist[i] = mean(samples$mz[valid])

    for (j in 1:length(valid)){ # Actualize intensity matrix
        I_matrix[i,samples$File_num[valid[j]]]= I_matrix[i,samples$File_num[valid[j]]]+samples$I[valid[j]]
      }}

  ########################
  ### Output new spectrum:
  ########################

  new_intensity = rowMeans(I_matrix)
  new_spectrum = cbind(avg_mzlist,new_intensity)
  colnames(new_spectrum)=NULL

  if (clean){ # Peak omni-present mass peaks
    kept = which(apply(I_matrix,1,function(x) all(x>0)))
    if (length(kept)>0){
      new_spectrum = new_spectrum[kept,]
      new_spectrum = matrix(new_spectrum,ncol=2)}
    }

  #write.table(I_matrix,"I_matrix.txt",col.names=F,row.names=F,sep="\t",dec=".")
  new_spectrum[,2]=new_spectrum[,2]/max(new_spectrum[,2])*100 # Normalized again to be safe

  return(new_spectrum)
}

############################
### Internal functions:
###########################

cut_mz_list<-function(mzlist,ppm_window){

  # The function cuts the ordered mass list according to user-defined tolerance window
  # It returns a feature list indicating which group each mass belongs to

  N=length(mzlist)

  # First cut: Neighbouring masses < ppm_window, temporary features

  f=1
  mz_feature=c(1,rep(0,N-1)) # which mass group is each mass

  for (k in 1:(N-1)){
    diff=(mzlist[k+1]-mzlist[k])/mzlist[k]*1000000
    if (diff<=ppm_window){
      mz_feature[k+1]=f}
    else{
      f=f+1
      mz_feature[k+1]=f}}

  # Second cut: Hierarchical clustering iniside each mass group

  feature_list=unique(mz_feature)
  new_mz_feature=rep("0",N) # Feature for second cut

  for (ff in feature_list){
    valid=which(mz_feature==ff)
    if (length(valid)>2){
      dis <- custom.dist(mzlist[valid], ppm_distance) # ppm distance
      discriminated<-cutree(hclust(dis,method="average"),h=ppm_window) # discriminated cluster of masses
      new_mz_feature[valid]=paste(ff,discriminated,sep="-")}

    else { # If only 1-2 features no HCT is provided
      new_mz_feature[valid]= paste(ff,"1",sep="-")}}

  return(new_mz_feature)
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

ppm_distance<-function(x,y){
  return(abs((x-y)/y*1000000))}

