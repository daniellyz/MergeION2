#' Spectral library query
#'
#' The function queries the input (reference) spectral library based on a query expression and/or (a) query MS/MS spectra or query library or LC-MS/MS file
#'
#' @param input_library Character or a list object. If character, file name with extension mgf, msp or RData. If list, must contain "complete", "consensus" and "network". 
#' @param query_ids Vector. Vectors of molecular ids that will be extracted from input library
#' @param query_expression Vector of characters or "". Vector of conditions used for querying the library. e.g. c("IONMODE=Positive","PEPMASS=325.19"). The left-hand side must match with the metadata items of the searched library.
#' @param query_spectrum Two-column data matrix. Two columns represent m/z and intensity of query tandem spectrum. At least 3 valid peaks should be provided. 
#' @param query_library Character or a list object. Used when query_spectrum = NULL. If character, file name with extension mgf, msp, Rdata, cdf, mz(x)ml. If list, must contain metadata and sp, or library_generator output ("complete", "consensus" and "network") 
#' @param params.search General parameters for searching spectral library
#' \itemize{
#'  \item{mz_search:}{ Numeric. Absolute mass tolerance in Da for fragment match.}
#'  \item{ppm_search:}{ Numeric. Mass tolerance in ppm for precursor mass match.} 
#'  \item{rt_search:}{ Numeric. Absolute retention time tolerance in second.}
#'  \item{rt_gap:}{ Numeric. Only used when querying a raw LC-MS/MS file. Retention time gap in second - when two scans both match with an input structure, they are both recorded as isomeric features of the same identifier if they are separated by a certain retention time gap. Please set it to 10000 if no isomeric feature is picked}
#' }
#' @param params.query.sp Parameters for matching an unknown spectra
#' \itemize{
#'   \item{prec_mz:}{ Numeric. Precursor mass of query spectrum.}
#'   \item{use_prec:}{ Boolean. If set to TRUE, precursor mass is used to "pre-query" the library.}
#'   \item{polarity:}{ Character. Ion mode of the query spectrum, either positive or negative.}
#'   \item{method:}{ Character. Similarity metrics.}
#'    \enumerate{
#'        \item{Precision:}{ Percentage of fragment or neutral loss matched in query spectrum.}
#'        \item{Recall:}{ Percentage of fragment or neutral loss matched in reference spectrum.}
#'        \item{F1:}{ Harmonic mean of precision and recall.}
#'        \item{Cosine:}{ COsine similarity score based on intensity vectors of fragments.}
#'        \item{Spearman:}{ Spearman similarity based on intensity ranks of fragments.}
#'        \item{MassBank:}{Similarity score used by MassBank using weighted pearson score}
#'        \item{NIST:}{Similarity score used by NIST using weighted pearson score}}
#'    \item{min_frag_match}{ Minimum matched peaks (or corresponding neutral losses) to allow a match.}
#'    \item{min_score}{ Minimum similarity score to be considered as a successful query.}
#'    \item{reaction_type:}{ Character. Either "Metabolic" and "Chemical". Type of transformation list used to annotate mass difference between connected features in molecular network.}
#' }
#'  
#' @return 
#' \itemize{
#'   \item{complete:}{ Complete spectra library filtered by both query experession and query spectrum}
#'   \item{consensus:}{ Consensus spectra library filtered by both query experession and query spectrum.}
#'   \item{network:}{ Molecular network filtered by both query experession and query spectrum. The query spectrum is also connected to the molecular network}
#'}
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#'
#' @examples
#'
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom plyr rbind.fill
#' 
#' @export
library_query<-function(input_library = NULL, query_ids = NULL, query_expression = "", query_spectrum = NULL, query_library = NULL,
            params.search = list(mz_search = 0.01, ppm_search = 10, rt_search = 15, rt_gap = 30), 
            params.query.sp = list(prec_mz = 10000, use_prec = T, polarity = "Positive", method = "Cosine", min_frag_match = 6, min_score = 0, reaction_type = "Metabolic")){
  
    options(stringsAsFactors = FALSE)
    options(warn=-1)
    gc()
    complete_selected = NULL
    consensus_selected = NULL
    network_selected = NULL
  
    mz_search = params.search$mz_search
    ppm_search = params.search$ppm_search
    rt_search = params.search$rt_search
    rt_gap  = params.search$rt_gap
  
    query_prec_mz  = params.query.sp$prec_mz
    query_use_prec = params.query.sp$use_prec
    query_polarity = params.query.sp$polarity
    query_method = params.query.sp$method
    query_min_frag = params.query.sp$min_frag_match
    query_min_score = params.query.sp$min_score
    if (is.null(query_min_score)){query_min_score = 0}
    query_reaction = params.query.sp$reaction_type
  
    if (query_reaction=="Metabolic"){
      data(reactionBio)
      #reactionList = read.csv("https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/reactionBio.txt", sep = "\t")
    }
    if (query_reaction=="Chemical"){
      data(reactionChem)
     # reactionList = read.csv("https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/reactionChem.txt", sep = "\t")
    }  
    
  ######################
  ### Filter Library ###
  ######################
    
  input_library = library_reader(input_library)
  
  # Filter on query ids:
    
  if (!is.null(query_ids)){
    if (!is.null(input_library$complete)){
        valid = which(input_library$complete$metadata$ID %in%query_ids)
        input_library$complete$metadata = input_library$complete$metadata[valid,,drop=FALSE]
        input_library$complete$sp = input_library$complete$sp[valid]
    }  
    if (!is.null(input_library$consensus)){
      valid = which(input_library$consensus$metadata$ID %in%query_ids)
      input_library$consensus$metadata = input_library$consensus$metadata[valid,,drop=FALSE]
      input_library$consensus$sp = input_library$consensus$sp[valid]
    }  
  }
    
  # Filter on query expressions:
    
    if (query_expression!=""){
      if (!is.null(input_library$complete$metadata)){
        input_library$complete = process_query(input_library$complete, query = query_expression, ppm_search, rt_search)$SELECTED
      }   
      if (!is.null(input_library$consensus$metadata)){
        input_library$consensus = process_query(input_library$consensus, query = query_expression, ppm_search, rt_search)$SELECTED
      }   
  }
    
  ###################################
  ### Load query sp or query file ###
  ###################################
  
  qs_library = NULL
  qs_network = NULL
    
  score_summary = c()
    
  if (!is.null(query_spectrum)){
    
     qs_library = list()
     
     tmp_qs_metadata = c("Spectrum", query_prec_mz, 0, "N/A", 2, nrow(query_spectrum))
     tmp_qs_metadata = as.data.frame(matrix(tmp_qs_metadata, nrow=1))
     colnames(tmp_qs_metadata) = c("ID", "PEPMASS", "RT", "FILENAME", "MSLEVEL", "TIC")
     qs_library$metadata = tmp_qs_metadata
     
     qs_library$sp = list()
     qs_library$sp[[1]] = query_spectrum
  }
  
  if (!is.null(query_library) & is.null(query_spectrum)){
    
    if (is.character(query_library)){
    
      if (tolower(file_ext(query_file)) %in% c("msp", "mgf")){
        tmp_qs = library_reader(query_file)
        qs_library = tmp_qs$complete
      }
    
      if (tolower(file_ext(query_file)) %in% c("cdf", "mzml", "mzxml")){
      # Start Automated feature detection
        tmp_qs = library_generator(input_library = NULL, lcms_files = query_file, metadata_file = NULL, 
            polarity = query_polarity, mslevel = 2, add.adduct = FALSE, processing.algorithm = "Default",
            params.search = list(mz_search = mz_search, ppm_search = ppm_search, rt_search = rt_search, rt_gap = rt_gap), 
            params.ms.preprocessing = list(normalized = TRUE, baseline = 1000, relative = 0.1, max_peaks = 500, recalibration = 0),
            params.consensus = list(consensus = TRUE, consensus_method = "consensus", consensus_window = 0.02),
            params.user = list(sample_type = "", user_name = "", comments = ""))
        qs_library = tmp_qs$complete
      }}
    
    if (is.list(query_library)){
      tmp_qs = library_reader(query_library)
      if (!is.null(tmp_qs$consensus)) {
        qs_library = tmp_qs$consensus
      } else {qs_library = tmp_qs$complete}
      if (!is.null(tmp_qs$network)) {qs_network = tmp_qs$network}
    }
  }
    
  if (!is.null(qs_library)){
    NQS = length(qs_library$sp) # Number of query spectra
    qs_metadata = qs_library$metadata
    qs_metadata$ID = paste0("Query_", qs_metadata$ID)
    qs_sp = qs_library$sp
    tmp_qs = list(complete = qs_library, consensus = qs_library, network = qs_network)
  } else {
    NQS = 0
    qs_metadata = NULL
    qs_sp = NULL
    tmp_qs = NULL
  }

  ####################
  ###Query spectrum ##
  ####################
  
  if (!is.null(input_library$complete)){ 
    complete_selected = input_library$complete
  }
  
  if (!is.null(input_library$consensus$metadata)){ 
    consensus_selected = input_library$consensus
    consensus_selected = process_query(consensus_selected, query = "MSLEVEL=2")$SELECTED
  } 
      
  if (!is.null(complete_selected) & is.null(consensus_selected) & !is.null(qs_sp)){
    
      message("Generating consensus MS/MS spectral library...")
       
      tmp_library = complete_selected
      
      input_library <-library_generator(input_library = tmp_library, mslevel = 2,
            params.search = list(mz_search = mz_search, ppm_search = ppm_search, rt_search = rt_search, rt_gap = 30), 
            params.ms.preprocessing = list(normalized = TRUE, baseline = 0, relative = 0, max_peaks = 200, recalibration = 0),
            params.consensus = list(consensus = TRUE, consensus_method = "consensus", consensus_window = mz_search*2),
            params.network = list(network = FALSE, similarity_method = query_method, min_frag_match = query_min_frag, min_score = 0.6, topK = 10, reaction_type = query_reaction, use_reaction = FALSE))
      consensus_selected =  input_library$consensus
  } 
      
  if (!is.null(consensus_selected) & !is.null(qs_sp)){
      
      id_selected = consensus_selected$metadata$ID
        
      tmp_library = list(complete = complete_selected, consensus = consensus_selected, network = input_library$network)
      
      for (jjj in 1:NQS){
        
        qs_mz = as.numeric(qs_metadata$PEPMASS[jjj])
        
        tmp_scores = process_similarity(qs_sp[[jjj]], 
              polarity = query_polarity, prec_mz = qs_mz, 
              use.prec = query_use_prec, input_library = tmp_library,  
              method = query_method, prec_ppm_search = ppm_search, 
              frag_mz_search = mz_search, min_frag_match = query_min_frag)
        tmp_scores = tmp_scores[tmp_scores$SCORES>=query_min_score,,drop=FALSE]
        
        if (!is.null(tmp_scores)){
          if (nrow(tmp_scores)>0){
            valid = which(round(tmp_scores$SCORES, 1) == max(round(tmp_scores$SCORES,1)))
            tmp_scores = tmp_scores[valid,,drop=FALSE]            
            tmp_scores = tmp_scores[which(!duplicated(tmp_scores$ID)),,drop=FALSE]
            idx = match(tmp_scores$ID, id_selected)
            tmp_ref = consensus_selected$metadata[idx,,drop=FALSE] 
            
            mdiff = round(abs(qs_mz - as.numeric(tmp_ref$PEPMASS)),3)
          #  tmp_annotation[i] = paste(tmp$ID[valid], collapse = ":")
           # tmp_mdiff[i] = paste(mdiff, collapse = ":")
          
            tmp_scores$QS = qs_metadata$ID[jjj]
            tmp_scores$MDIFF = mdiff
            
            score_summary = rbind.data.frame(score_summary, tmp_scores)
        }}
    }
      
    if (!is.null(score_summary)){
        if (nrow(score_summary)>0){
          
          # Update 
          
          score_ids = unique(score_summary$ID)
          query_ids = unique(score_summary$QS)
          
          idx = match(score_ids, complete_selected$metadata$ID)
          complete_selected = list(metadata =  complete_selected$metadata[idx,,drop=FALSE], sp = complete_selected$sp[idx])
          
          idx = match(score_ids, consensus_selected$metadata$ID)
          consensus_selected = list(metadata =  consensus_selected$metadata[idx,,drop=FALSE], sp = consensus_selected$sp[idx])
          
      }}  
  }
  
  ##########################################
  ### Build Query-Library Hybrid Network ###
  ##########################################
  
  if (!is.null(input_library$network$network) & !is.null(score_summary)){
    
    # Update profile matrix::
    
    db_profile = input_library$network$db_profile[,idx,drop=FALSE]
    valid = which(apply(db_profile, 1, function(x) sum(x>0))>0) # Remove all zero features
    db_profile = db_profile[valid,,drop=FALSE]
    db_feature = input_library$network$db_feature[valid,,drop=FALSE]
  
    # Update nodes:
    
    idx = match(score_ids, input_library$network$nodes$ID)
    lib_nodes = input_library$network$nodes[idx,,drop=FALSE] # Nodes of library matched

    query_nodes = rbind.fill(lib_nodes, qs_metadata)
    
    # Update network:

    # Network of library matched

    idx2a = match(score_ids, input_library$network$network[,1])
    idx2b = match(score_ids, input_library$network$network[,2])
    
    idx2  = intersect(idx2a, idx2b)
    idx2 = idx2[!is.na(idx2)]

    lib_network = input_library$network$network[idx2,,drop=FALSE]
    lib_network = lib_network[,1:4]

    # Network of query spectra
    
    test_network = tmp_qs$network$network
    
    if (NQS>1 & is.null(test_network)){
    
      test_network = process_lib2network(tmp_qs, networking = T, polarity = query_polarity,
                params.search = list(mz_search = mz_search, ppm_search = ppm_search),
                params.similarity = list(method = query_method, min.frag.match = query_min_frag, min.score = query_min_score),
                params.network = list(topK = 10, max.comp.size = 50, reaction.type = query_reaction, use.reaction = F))
    }
    
    if (!is.null(test_network)){
      test_network = test_network[,1:4,drop=FALSE]
      test_network[,1] = paste0("Query_", test_network[,1])
      test_network[,2] = paste0("Query_", test_network[,2])
    }
    
    # Network of query-library matches
    
    qs_network = cbind.data.frame(ID1 = score_summary$QS, ID2 = score_summary$ID, 
              MS2.Matches = score_summary$PEAK.MATCHES, MS2.Similarity = score_summary$SCORES) 

    # Combine three networks
    
    query_network = rbind.fill(lib_network, test_network, qs_network)
    
    # Extra info to network
    
    NE = nrow(query_network)
    reaction_annotated = rep("N/A", NE)
    reaction_formula = rep("N/A", NE)
    node_id = as.character(query_nodes$ID) 

    for (k in 1:NE){
           
        II1 = which(node_id == as.character(query_network$ID1[k]))[1]
        MZ1 = as.numeric(query_nodes$PEPMASS[II1])
      
        II2 = which(node_id == as.character(query_network$ID2[k]))[1]
        MZ2 = as.numeric(query_nodes$PEPMASS[II2])
           
        MDiff = abs(MZ1 - MZ2)
        MDiff_error = ppm_distance1(MDiff, reactionList$Mdiff)

        if (min(MDiff_error)<=ppm_search){
          
          ind = which.min(MDiff_error)[1]
          reaction_annotated[k] = reactionList$Reaction.Name[ind]
          reaction_formula[k] = reactionList$Formula[ind]
        }
    }
           
    query_network = cbind.data.frame(query_network, reaction = reaction_annotated, reaction_formula = reaction_formula)
    network_selected = list(db_profile = db_profile, db_feature = db_feature, nodes = query_nodes, network = query_network)
  }
  
  #############
  ###Output ###
  #############
  
  if (NQS>1){
    tmp_annotation = rep("0", NQS)
    tmp_scores = rep(0, NQS)
    tmp_mdiff = rep("0", NQS)

    for (i in 1:NQS){
      
      tmp_id = paste0("Query_", qs_library$metadata$ID[i])
      idx = which(score_summary$QS == tmp_id)
      if (length(idx)>0){
        sss = score_summary[idx,]
        tmp_annotation[i] = paste(sss$ID, collapse = ":")
        tmp_scores[i] = paste(round(sss$SCORES,2), collapse = ":")
        tmp_mdiff[i] = paste(round(sss$MDIFF,3), collapse = ":")
      }
    }
    
    qs_library$metadata$ANNOTATION_ANALOGUE = tmp_annotation
    qs_library$metadata$SCORE_ANALOGUE = tmp_scores
    qs_library$metadata$MDIFF_ANALOGUE = tmp_mdiff
    qs_library$metadata$PEPMASS = round(as.numeric(qs_library$metadata$PEPMASS), 3)
    output_library = list(complete = qs_library, consensus = qs_library, network = network_selected)
    return(output_library)
  }
    
  if (NQS==1){
    consensus_selected$metadata$SCORE_MERGEION = paste(score_summary$SCORES, collapse = ":")
    output_library = list(complete = complete_selected, consensus = consensus_selected, network = network_selected)
    return(output_library)
  }
  
  if (NQS==0 & !is.null(complete_selected)){
    output_library = list(complete = complete_selected, consensus = consensus_selected, network = network_selected)
    return(output_library)
  }
  
}

##########################
### Internal functions:###
##########################

ppm_distance1<-function(x,y){
  x = as.numeric(x)
  y = as.numeric(y)
  if (y>100){
    ppm = abs((x-y))/y*1000000
  } else {
    ppm = abs(x-y)
    ppm[ppm<=0.01]=0
    ppm[ppm>0.01]=50
  }
  return(ppm)
}
