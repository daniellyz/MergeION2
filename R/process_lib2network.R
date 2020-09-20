#' Create feature-based molecular networks from a conesnsus spectra library
#'
#' Function used by library_generator to create molecular networks
#'
#' @param input_library A list object. Must contain consensus library.
#' @param networking Boolean. TRUE to go through the entire molecular networking process. FALSE if only spectra alignment is performed, necessary for spectral library searching
#' @param polarity character. Either "Positive" or "Negative". Ion mode of the LC-MS/MS file. 
#' @param params.screening Parameters for feature screening and MS2 spectra pre-processing from raw_data_file:
#' \itemize{
#'  \item{baseline:}{ Numeric. Absolute intensity threshold that is considered as a mass peak and written into the library.}
#'  \item{relative:}{ Numeric between 0 and 100. Relative intensity threshold of the highest peak in each spectrum, peaks above both absolute and relative thresholds are saved in the library.}
#'  \item{max_peaks:}{ Integer higher than 3. Maximum number of peaks kept per spectrum from the highest peak.}
#' }
#' @param params.search List of parameters for feature screening, combining precursor ions and fragments in the input file, as well as for searching in the spectral library. The list must contain following elements:
#' \itemize{
#'  \item{mz_search:}{ Numeric. Absolute mass tolerance in Da.}
#'  \item{ppm_search:}{ Numeric. Absolute mass tolerance in ppm.} 
#' }
#' @param params.similarity Parameters for MS/MS spectral similarity determination, used for both molecular networking and spectral library search.
#' \itemize{
#'  \item{method:}{ Characeter.Similarity metrics for networking and spectral library search. Must be "Matches", "Dot", "Cosine", "Spearman", "MassBank", "NIST". Please check function library_query for more details.}
#'  \item{min.frag.match:}{ Integer. Minimum number of common fragment ions (or neutral losses) that are shared to be considered for spectral similarity evaluation. We suggest setting this value to at least 6 for statistical meaningfulness.
#'  \item{min.score:}{ Numeric between 0 and 1. Minimum similarity score to annotate an unknown feature with spectral library or to connect two unknown features because they are similar. It does NOT affect method = "Matches".}
#' }
#'@param params.network Parameters for post-filtering and annotation of network edges: based on feature correlation (if feature quantification table is provided) and mass difference
#' \itemize{
#'  \item{topK:}{ Integer higher than 0. For networking, the edge between two nodes are kept only if both nodes are within each other's TopK most similar nodes. For example, if this value is set at 20, then a single node may be connected to up to 20 other nodes. Keeping this value low makes very large networks (many nodes) much easier to visualize. We suggest keeping this value at 10.}
#'  \item{reaction.type:}{ Character. Either "Metabolic" and "Chemical". Type of transformation list used to annotate mass difference between connected features in molecular network.}
#'  \item{use.reaction:}{ Boolean. TRUE if keep only edges whose mass difference can be annotated to known metabolic or chemical reactions.}
#' }
#'  
#' @export
#'
process_lib2network<-function(input_library, networking = T, polarity = c("Positive", "Negative"),
                  params.screening = list(baseline = 1000, relative = 0.01, max_peaks = 200),
                  params.search = list(mz_search = 0.005, ppm_search = 10),
                  params.similarity = list(method = "Cosine", min.frag.match = 6, min.score = 0.6),
                  params.network = list(topK = 10, reaction.type = "Metabolic", use.reaction = TRUE)){
  
  options(stringsAsFactors = FALSE)
  options(warn=-1)

  ####################
  ###  Input Check ###
  ####################

  input_library = library_reader(input_library)
  
  complete_library = input_library$complete
  consensus_library0 = input_library$consensus

  if (is.null(consensus_library0)){stop("No consensus spectra available!")}
  consensus_library = process_query(consensus_library0, query = "MSLEVEL=2")$SELECTED
  if (is.null(consensus_library)){stop("No MS2 in consensus spectra available!")}
  
  NI = nrow(consensus_library$metadata)
  
  metadata = consensus_library$metadata
  splist = consensus_library$sp
  IDList = metadata$ID
  MZList = as.numeric(metadata$PEPMASS)
  
  baseline = params.screening$baseline
  relative = params.screening$relative
  max_peaks = params.screening$max_peaks 
  
  mz_search = params.search$mz_search
  ppm_search = params.search$ppm_search

  sim.method = params.similarity$method
  min.frag.match =  params.similarity$min.frag.match
  min.score = params.similarity$min.score
  
  topK = params.network$topK
  reaction.type = params.network$reaction.type
  use.reaction = params.network$use.reaction

  if (!(reaction.type %in% c("Chemical", "Metabolic"))){
    stop("Reaction type for annotating network mass differences must be Chemical or Metabolic!")
  }
  
  ###################################
  ### Transform spectra to matrix ###
  ###################################

  library_matrix = matrix_generator(consensus_library, mz_window = mz_search)
  
  new_nodes = metadata
  new_network = c()
  
  #######################################
  ### Spectral similarity calculation ###
  #######################################
  
  if (networking){
  
    message("Generating molecular network...")
    
    for (i in 1:(NI-1)){
  
      temp_spectrum = splist[[i]]
      temp_library = list(complete = NULL, consensus = NULL, network = NULL)
      
      # Search part of library
      
      lib_range = (i+1):NI
      temp_library$consensus$metadata = consensus_library$metadata[lib_range,,drop=FALSE]
      temp_library$consensus$sp = consensus_library$sp[lib_range]
      temp_library$network$db_profile =  library_matrix$db_profile[,lib_range,drop=FALSE] 
      temp_library$network$db_feature =  library_matrix$db_feature
      
      temp_scores = process_similarity(query_spectrum = temp_spectrum, polarity = polarity, prec_mz = MZList[i], use.prec = FALSE, input_library = temp_library,
            method = sim.method, prec_ppm_search = ppm_search, frag_mz_search = mz_search, min_frag_match = min.frag.match)

      if (!is.null(temp_scores)){
        NSC = nrow(temp_scores)
        temp_network = cbind(ID1 = rep(IDList[i],NSC), ID2 = temp_scores[,1], MS2.Similarity = round(temp_scores[,2],3))  
        new_network = rbind.data.frame(new_network, temp_network)
      }
    }

    #######################
    ### Top K filtering ###
    #######################
  
    if (!is.null(new_network)){
      new_network = mutual_filter(new_network, topK = topK)
      if (nrow(new_network)==0){new_network = NULL}
    }
    
    ##################################
    ### Mass difference annotation ###
    ##################################
  
    if (!is.null(new_network)){
  
      if (reaction.type=="Metabolic"){
        reactionList = read.csv("https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/reactionBio.txt", sep = "\t")
      }
      if (reaction.type=="Chemical"){
        reactionList = read.csv("https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/reactionChem.txt", sep = "\t")
      }
    
      reaction_annotated = rep("N/A", nrow(new_network))
      reaction_formula = rep("N/A", nrow(new_network))
  
      for (k in 1:nrow(new_network)){
        II1 = which(as.character(metadata$ID) == as.character(new_network$ID1[k]))[1]
        II2 = which(as.character(metadata$ID) == as.character(new_network$ID2[k]))[1]
        MZ1 = as.numeric(metadata$PEPMASS[II1])
        MZ2 = as.numeric(metadata$PEPMASS[II2])
      
        MDiff = abs(MZ1 - MZ2)
        MDiff_error = ppm_distance(MDiff, reactionList$Mdiff)
        
        if (min(MDiff_error)<=ppm_search){
          ind = which.min(MDiff_error)[1]
          reaction_annotated[k] = reactionList$Reaction.Name[ind]
          reaction_formula[k] = reactionList$Formula[ind]
        }
    }
    
    new_network = cbind.data.frame(new_network, reaction = reaction_annotated, reaction_formula = reaction_formula)
  }
  
  ######################################
  ### Final filtering and outputing  ###
  ######################################
  
  if (!is.null(new_network)){
    
      new_network = new_network[new_network[,3]>=min.score,,drop=FALSE]
    
      if (use.reaction){
        new_network = new_network[which(new_network$reaction!="N/A"),,drop=FALSE]
      }
      
      if (nrow(new_network)==0){new_network = NULL}
    }
  }
  
  output_network = list(db_profile = library_matrix$db_profile, 
                        db_feature = library_matrix$db_feature, 
                        nodes = new_nodes, network = new_network)
  
  output_library = list(complete = complete_library, consensus = consensus_library0, network = output_network)
  
  return(output_library)
}
                  
                  
############################################
### Transforming input library to matrix:###
############################################

matrix_generator<-function(input_library, mz_window = 0.02){
  
  # Generate profile matrix (N identical structures and M features) and feature list from input library
  
  # Some exceptions:
  
  if (is.null(input_library$metadata)){return(NULL)}
  
  if (nrow(input_library$metadata)==0){return(NULL)}
  
  if (nrow(input_library$metadata)==1){
    
    frags =  input_library$sp[[1]][,1]
    ints =  input_library$sp[[1]][,2]
    nls = as.numeric(input_library$metadata$PEPMASS)-frags
    FID = paste0("Frag_", 1:length(frags))
    NID = paste0("Nloss_", 1:length(nls))
    
    sp_profile =  nl_profile = t(t(ints))
    colnames(sp_profile) = colnames(nl_profile) = input_library$metadata$ID
    sp_feature =  cbind.data.frame(FID = FID, Mass = frags)
    nl_feature =  cbind.data.frame(FID = NID, Mass = nls)
    colnames(sp_feature) = colnames(nl_feature) = c("ID", "Mass")
    
    db_profile = rbind(sp_profile, nl_profile)
    db_feature = cbind(rbind(sp_feature, nl_feature), Type = c(rep("Frag", nrow(sp_feature)), rep("Nloss", nrow(nl_feature))))
    
    return(list(db_profile = db_profile, db_feature = db_feature))
  }
  
  ### Otherwise:
  
  splist = input_library$sp
  metadata = input_library$metadata
  IDlist = metadata$ID
  NM = nrow(metadata)
  
  ### Neutral loss spectra:
  
  nllist = list()
  for (i in 1:NM){
    sp = splist[[i]]
    prec.mz = as.numeric(metadata$PEPMASS[i])
    nl = cbind(prec.mz - sp[,1], sp[,2])
    nl = nl[nl[,1]>2,,drop=FALSE]
    nl = nl[order(nl[,1]),,drop=FALSE]
    nllist[[i]] = nl
  }  
  
  ### Align spectra

  sp_aligned = average_spectrum(splist, mz_window = mz_window)
  sp_profile = sp_aligned$I_matrix
  FID = paste0("Frag_", 1:nrow(sp_profile))
  rownames(sp_profile) = FID
  colnames(sp_profile) = IDlist
  sp_feature = data.frame(FID = FID, Mass = sp_aligned$new_spectrum[,1])
  
  nl_aligned = average_spectrum(nllist, mz_window = mz_window)
  nl_profile = nl_aligned$I_matrix
  NID = paste0("Nloss_", 1:nrow(nl_profile))
  rownames(nl_profile) = NID
  colnames(nl_profile) = IDlist
  nl_feature = data.frame(NID = NID, Mass = nl_aligned$new_spectrum[,1])
  
  ### Reduce Existing Library

  colnames(sp_feature) = colnames(nl_feature) = c("ID", "Mass")
  
  db_profile = rbind(sp_profile, nl_profile)
  db_feature = cbind(rbind(sp_feature, nl_feature), Type = c(rep("Frag", nrow(sp_feature)), rep("Nloss", nrow(nl_feature))))
  
  valid = which(apply(db_profile, 1, sum)>0)
  db_profile = db_profile[valid,,drop=FALSE]
  db_feature = db_feature[valid,,drop=FALSE]
  
  return(list(db_profile = db_profile, db_feature = db_feature))
}  

######################
### Filter network:###
######################

mutual_filter <- function(network, topK = 10){
  
  # The function keeps edges in the network if both nodes are among the topK of each other
  
  NR = nrow(network)
  temp_id = unique(c(network[,1], network[,2]))
  NI = length(temp_id)
  
  # Transformation to matrix:
  
  MMM = matrix(0, NI, NI)  
  colnames(MMM) = rownames(MMM) = temp_id
  
  from_ind = match(network[,1], temp_id)
  to_ind = match(network[,2], temp_id) # Index in the temp_ind
  
  for(x in 1:NR){
    MMM[from_ind[x], to_ind[x]] = as.numeric(network[x,3])
    MMM[to_ind[x], from_ind[x]] = as.numeric(network[x,3])
  }
  
  # Take topK of each row
  
  top_value_row = apply(MMM, 1, function(x) min(head(sort(x, decreasing=T), topK))) # Lowest value of topK
  top_value_column = apply(MMM, 2, function(x) min(head(sort(x, decreasing=T), topK)))
  
  for (x in 1:NI){
    tmp = MMM[x,]
    tmp[tmp<top_value_row] = 0
    MMM[x,] = tmp
    tmp = MMM[,x]
    tmp[tmp<top_value_column] = 0
    MMM[,x] = tmp
  }  
  
  MMM[lower.tri(MMM)] = 0 # Set to 0 due to symmetry
  
  # Transform matrix back to network
  
  valid = which(MMM>0, arr.ind = T)
  from_id = temp_id[valid[,1]]
  to_id = temp_id[valid[,2]]
  sim_score = MMM[valid]
  
  network_filtered = cbind.data.frame(ID1 = from_id, ID2 = to_id, MS2.Similarity = MMM[valid])
  
  return(network_filtered)
}

##########################
### Internal functions:###
##########################

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