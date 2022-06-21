# Project: MergeION2
# 
# Author: yzhan482
###############################################################################

# rex from rex package
# dplyr filter %<>% %>% group_by

interactive_visulizer <- function(input_library, id, 
		type = c("complete", "consensus"),
		query_spectrum=NULL, collpase = TRUE){
	
stopifnot(!missing(input_library), !missing(id))	
	
input_library = library_reader(input_library)
type = match.arg(type)


if (!type %in% c("complete", "consensus")){
	stop("Please choose the visualized library type between complete and consensus!")
}

if (is.null(input_library$consensus)){

	message("No consensus library is available! The entire library is used for visualization!")
	type = "complete"
}


input_library <- switch(type, 
		"complete" = input_library$complete,
		"consensus" = input_library$consensus
		)



if (nrow(metadata)==0){stop("The library does not contain the selected ID!")}

## query based on ID
#library1 = process_query(input_library, query = paste0("ID =", id))$SELECTED


# query on multiple ID
	
idsQuery <- paste0("ID@xxx@", 
		paste0(sapply(id, function(x){paste0("(^",rex(x),"$)")}), collapse = "|")
)


library1 = process_query(library0 = input_library, query = idsQuery)$SELECTED

metadata = library1$metadata
spectrum_list = library1$sp

if (nrow(metadata)==0) stop("The library does not contain the selected ID(s)!")

# collapse all spectrum into one data frame
spectrum <- do.call("rbind", spectrum_list) %>% as_tibble()
colnames(spectrum) <- c("m/z", "Intensity")

# prepare all metadata
metadataExpand <- metadata[rep(1:nrow(metadata), times= unlist(lapply(spectrum_list, nrow))), ]
#dim(metadataExpand); dim(spectrum)

# combine spectrum and meta data	
plotData <- bind_cols(spectrum, metadataExpand)

# for each ID, MSLELVE select the largest scan
plotData %<>% 
		#.$SCANS %>% unique
		group_by(ID, MSLEVEL) %>% filter(SCANS == max(SCANS)) %>%
		filter(Intensity > max(Intensity)*0.05) %>%
		slice(1:100) %>%
		ungroup() 

p <- ggplot(data = plotData, aes(x = `m/z`, y = Intensity, label = format(Intensity, digits = 4), color = ID))
p <- p + geom_text(vjust =  0.1) +
		geom_segment(aes(xend =  `m/z`, yend = 0), size = 1, lineend = "butt") +
		facet_grid(MSLEVEL~., scale = "free_x", labeller = labeller(.rows = label_both)) +
		theme_bw() 
		
library(plotly)

ggplotly(p)
				
}