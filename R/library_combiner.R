#' Combine multiple libraries
#'
#' The function generates a combined library from several individual libraries 
#'
#' @param library1 First library. A list of object from library_generator, or library file names (mgf, rdata or msp format)
#' @param library2 Second library 
#' 
#' @importFrom plyr rbind.fill
#' 
#' @export
#'
#' @examples
#'
#' \dontrun{data(DRUG_THERMO_LIBRARY)
#' combined_lib = library_combiner(library1, library2)}
#'
#'
library_combiner<-function(...){
  
    all_libraries = list(...)
    NL = length(all_libraries)

    combined_metadata1 = c()
    combined_sp1 = c()
    
    combined_metadata2 = c()
    combined_sp2 = c()
    
    for (i in 1:NL){
      
      tmp_library = library_reader(all_libraries[[i]])

      # Complete:
      
      combined_metadata1 = rbind.fill(combined_metadata1, tmp_library$complete$metadata)
      combined_sp1 = c(combined_sp1, tmp_library$complete$sp)
      
      # Consensus:

      if (!is.null(tmp_library$consensus$metadata)){
        combined_metadata2 = rbind.fill(combined_metadata2, tmp_library$consensus$metadata)
        combined_sp2 = c(combined_sp2, tmp_library$consensus$sp)
      }
    }
    
    combined_library = list(complete = list(metadata = combined_metadata1, sp = combined_sp1),
                            consensus = list(metadata = combined_metadata2, sp = combined_sp2),
                            network = NULL)
    
    
    return(combined_library)
}
