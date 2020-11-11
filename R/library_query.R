#' Spectral library query
#'
#' The function queries the input spectral library based on a query expression and/or (a) query MS/MS spectra or raw LC-MS/MS file
#'
#' @param input_library Character or a list object. If character, file name with extension mgf, msp or RData. If list, must contain "complete", "consensus" and "network". 
#' @param query_expression Vector of characters. Vector of conditions used for querying the library. e.g. c("IONMODE=Positive","PEPMASS=325.19"). The left-hand side must match with the metadata items of the searched library.
#' @param query_spectrum  Two-column data matrix. Two columns represent m/z and intensity of query tandem spectrum. At least 3 valid peaks should be provided. 
#' @param query_file Character. The file name of query spectra collection, must in mgf, msp or mz(X)ML, cdf format. 
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
#'        \item{Matches:}{ Conting number of fragment or neutral loss matches}
#'        \item{Dot:}{ Dot similarity score based on intensity vectors of fragments.}
#'        \item{Cosine:}{ COsine similarity score based on intensity vectors of fragments.}
#'        \item{Spearman:}{ Spearman similarity based on intensity ranks of fragments.}
#'        \item{MassBank:}{Similarity score used by MassBank using weighted pearson score}
#'        \item{NIST:}{Similarity score used by NIST using weighted pearson score}}
#'    \item{min_frag_match}{ Minimum matched peaks (or corresponding neutral losses) to allow a match.}
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
#'   
#' @export

library_query<- function(input_library = NULL, query_expression = "IONMODE = Positive", query_spectrum = NULL, query_file = NULL,
                         params.search = list(mz_search = 0.01, ppm_search = 10, rt_seach = 15, rt_gap = 30), 
                         params.query.sp = list(prec_mz = 100, use_prec = T, polarity = "Positive", method = "Cosine", min_frag_match = 6, reaction_type = "Metabolic")){
  
  options(stringsAsFactors = FALSE)
  options(warn=-1)
  gc()
  complete_selected = NULL
  consensus_selected = NULL
  network_selected = NULL
  
  ############
  ### Read ###
  ############
  
  input_library = library_reader(input_library) 
  
  qs_list = list()
  
  if (!is.null(query_spectrum)){qs_list[[1]] = query_spectrum}
  
  if (!is.null(query_file)){
    if (tolower(file_ext(query_file)) %in% c("msp", "mgf")){
      tmp_qs = library_reader(query_file)
      qs_list = tmp_qs$library_complete
    }
  }
  
  mz_search = params.search$mz_search
  ppm_search = params.search$ppm_search
  rt_search = params.search$rt_seach
  rt_gap  = params.search$rt_gap
  
  query_prec_mz  = params.query.sp$prec_mz
  query_use_prec = params.query.sp$use_prec
  query_polarity = params.query.sp$polarity
  query_method = params.query.sp$method
  query_min_frag = params.query.sp$min_frag_match
  query_reaction = params.query.sp$reaction_type
    
  if (query_reaction=="Metabolic"){
    reactionList = read.csv("https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/reactionBio.txt", sep = "\t")
  }
  if (query_reaction=="Chemical"){
    reactionList = read.csv("https://raw.githubusercontent.com/daniellyz/MergeION2/master/inst/reactionChem.txt", sep = "\t")
  }
  
  ################################
  ###Query 1: Complete Library ###
  ################################
    
  if (!is.null(input_library$complete)){
      if (query_expression!=""){
        complete_selected = process_query(input_library$complete, query = query_expression, ppm_search, rt_search)$SELECTED
      } else {
        complete_selected = input_library$complete
      }
  } 
  
  #################################
  ###Query 2: Consensus Library ###
  #################################
      
  if (!is.null(complete_selected) & is.null(input_library$consensus) & !is.null(qs_list)){
    
       message("Generating consensus MS/MS spectral library...")
       
       tmp_library = list(complete = complete_selected)
       
       input_library <-library_generator(input_library = tmp_library, mslevel = 2,
                 params.search = list(mz_search = mz_search, ppm_search = ppm_search, rt_seach = rt_search, rt_gap = 30), 
                 params.ms.preprocessing = list(normalized = TRUE, baseline = 1000, relative = 0.1, max_peaks = 200, recalibration = 0),
                 params.consensus = list(consensus = TRUE, consensus_method = "consensus", consensus_window = mz_search*2),
                 params.network = list(network = TRUE, similarity_method = query_method, min_frag_match = query_min_frag, min_score = 0.6, topK = 10, reaction_type = query_reaction, use_reaction = FALSE))
   } 
      
  if (!is.null(input_library$consensus)){
      
      consensus_ids = input_library$consensus$metadata$ID
      
      if (query_expression!=""){
        consensus_selected = process_query(input_library$consensus, query = query_expression, ppm_search, rt_search)$SELECTED
      } else {
        consensus_selected = input_library$consensus
      }

      id_selected = consensus_selected$metadata$ID
        
      if (length(query_spectrum)>0){
        
        tmp_library = list(complete = complete_selected, consensus = consensus_selected, network = input_library$network)
        
        tmp_scores = process_similarity(query_spectrum, polarity = query_polarity, prec_mz = query_prec_mz, use.prec = query_use_prec, input_library = tmp_library,  
                    method = query_method, prec_ppm_search = ppm_search, frag_mz_search = mz_search, min_frag_match = query_min_frag)
        
        id_selected = tmp_scores$ID

        idx = match(id_selected, complete_selected$metadata$ID)
        complete_selected = list(metadata =  complete_selected$metadata[idx,,drop=FALSE], sp = complete_selected$sp[idx])
        
        idx = match(id_selected, consensus_selected$metadata$ID)
        consensus_selected = list(metadata =  consensus_selected$metadata[idx,,drop=FALSE], sp = consensus_selected$sp[idx])
        consensus_selected$metadata$SCORE_MERGEION = tmp_scores$SCORES
    }
  }
  
  #######################
  ###Query 3: Network ###
  #######################
  
  if (!is.null(input_library$network$network) & length(id_selected)>0){
        
     idx = match(id_selected, consensus_ids)

     # Update profile matrix::
        
     db_profile = input_library$network$db_profile[,idx,drop=FALSE]
     valid = which(apply(db_profile, 1, function(x) sum(x>0))>0) # Remove all zero features
     db_profile = db_profile[valid,,drop=FALSE]
     db_feature = input_library$network$db_feature[valid,,drop=FALSE]
        
     # Update Nodes:
        
     idx1 = match(id_selected, input_library$network$nodes$ID)
     nodes = input_library$network$nodes[idx1,,drop=FALSE]
        
     if (!is.null(query_spectrum)){
         nodes$MERGEION_SIMILARITY_SCORE =  tmp_scores$SCORES
         query_node = nodes[1,,drop=FALSE]  
         query_node[1,] = "N/A"
         query_node$PEPMASS = query_prec_mz
         query_node$IONMODE = query_polarity
         query_node$ID = "Query_Spectrum"
         query_node$SCANS = -999
         nodes = rbind.data.frame(nodes, query_node)
      }
        
     # Update netwerk:
        
     idx2a = match(id_selected, input_library$network$network[,1])
     idx2b = match(id_selected, input_library$network$network[,2])
     idx2  = union(idx2a, idx2a)
     idx2 = idx2[!is.na(idx2)]
     network = input_library$network$network[idx2,,drop=FALSE]
        
     if (!is.null(query_spectrum)){
       
         query_network = cbind.data.frame(ID1 = "Query_Spectrum", ID2  = tmp_scores$ID, MS2.Similarity = tmp_scores$SCORES)
         NE = nrow(query_network)
         
         reaction_annotated = rep("N/A", NE)
         reaction_formula = rep("N/A", NE)
         
         node_id = as.character(nodes$ID) 
         node_id = node_id[-length(node_id)]
         
         for (k in 1:NE){
           
           II2 = which(node_id == as.character(query_network$ID2[k]))[1]
           MZ2 = as.numeric(nodes$PEPMASS[II2])
           
           MDiff = abs(query_prec_mz - MZ2)
           MDiff_error = ppm_distance(MDiff, reactionList$Mdiff)
           
           if (min(MDiff_error)<=ppm_search){
             ind = which.min(MDiff_error)[1]
             reaction_annotated[k] = reactionList$Reaction.Name[ind]
             reaction_formula[k] = reactionList$Formula[ind]
           }
         }
           
         query_network = cbind.data.frame(query_network, reaction = reaction_annotated, reaction_formula = reaction_formula)
         network = rbind.data.frame(network, query_network)
    }

    network_selected = list(db_profile = db_profile, db_feature = db_feature, nodes = nodes, network = network)
  }
  
  #############
  ###Output ###
  #############
     
  output_library = list(complete = complete_selected, consensus = consensus_selected, network = network_selected)
  
  return(output_library)
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