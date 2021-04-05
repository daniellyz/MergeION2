#' Searching and managing the spectral library
#' 
#' The function is used by library_query() to select scans based on a query expression.
#' 
#' @export
#'
#' @importFrom stringr str_replace_all fixed
#'
process_query<-function(library0, query = "", ppm_search = 20, rt_search = 12){

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  #####################################
  ### Reading from spectral library:###
  #####################################

  metadata = library0$metadata
  spectrum_list = library0$sp

  prec_mz = as.numeric(metadata$PEPMASS)
  prec_rt = as.numeric(metadata$RT)
  
  #############################################
  ### Run query expressions and find indexes ##
  #############################################

  if (!is.character(query)){
    stop("Query expression is not valid!")}

  if (query!=""){

    indexes_list = list()
    NI = 0

    for (eps in query){
      eps1 =  str_replace_all(eps,fixed(" "),"") # Remove white space

    # Search pepmass and rt:
      
      if (startsWith(eps1,"PEPMASS=")){
        target_mass = as.numeric(strsplit(eps1,"=")[[1]][2])
        if (!is.na(target_mass)){
          ppm_list = ppm_distance(prec_mz, target_mass)
          indexes = which(ppm_list<=ppm_search)
        }
      } else if (startsWith(eps1,"RT=")){
        target_rt = as.numeric(strsplit(eps1,"=")[[1]][2])
        if (!is.na(target_rt)){
          rtdev_list = abs(target_rt*60-prec_rt*60)
          indexes = which(rtdev_list<=rt_search)
      }} else {

    # Search other things:
      target_variable = strsplit(eps1,"=")[[1]][1]
      target_value = strsplit(eps1,"=")[[1]][2]
      cid = which(colnames(metadata) == target_variable)
      if (length(cid)==1){
        indexes = which(metadata[,cid]==target_value)}}

    # Add valid indexes:

      if (length(indexes)>0){
          NI = NI + 1
          indexes_list[[NI]] = indexes}
    }

    indexes_list = Reduce(intersect,indexes_list)
   }

  # Output results:

  NN = 1:length(spectrum_list)
  left_list = setdiff(NN, indexes_list)

  SELECTED_LIBRARY = LEFT_LIBRARY = library0
  SELECTED_LIBRARY$sp = library0$sp[indexes_list]
  SELECTED_LIBRARY$metadata = library0$metadata[indexes_list,,drop=FALSE]
  LEFT_LIBRARY$sp = library0$sp[left_list]
  LEFT_LIBRARY$metadata = library0$metadata[left_list,,drop=FALSE]

  return(list(SELECTED = SELECTED_LIBRARY,
              LEFT = LEFT_LIBRARY))
}

#########################
### Internal functions:##
#########################

# ppm error calculation:

ppm_distance<-function(x,y){
  x = as.numeric(x)
  y = as.numeric(y)
  if (y>100){
    ppm = abs((x-y))/y*1000000
  } else {
    ppm = abs(x-y)
    ppm[ppm<0.01]=0
  }
  return(ppm)
}

