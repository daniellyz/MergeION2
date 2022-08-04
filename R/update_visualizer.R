
# input_library <- library_reader("//10.157.249.83/share/mergeion/JANSSEN_COMPLETE.RData")
# input_library_sub <- library_reader("//10.157.249.83/share/mergeion/JANSSEN_CONSENSUS_POS_SUBSTRUCTURE.RData")
# id =  input_library$complete$metadata$ID[1:100]
# type = "complete"

#' Function to generate data for interactive visulization
#' @inheritParams library_visualizer
#' @importFrom rex rex
#' @importFrom dplyr fitler group_by bind_cols
#' @importFrom magrittr %>%  %<>%
#' @return a dataframe with m/z, intensity, ID and all the metadata columns from \code{input_library} metadata
#' @author Yingjie Zhang
#' @export 

dataForPlotly <- function(input_library, 
		id, 
		type = c("complete", "consensus"),
		query_spectrum=NULL){
	
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
	
	
	input_librarySub <- switch(type, 
			"complete" = input_library$complete,
			"consensus" = input_library$consensus
	)
	
	
	
# query on multiple ID
	
	idsQuery <- paste0("ID@xxx@", 
			paste0(sapply(id, function(x){paste0("(^",rex(x),"$)")}), collapse = "|")
	)
	
	
	library1 <- process_query(library0 = input_librarySub, query = idsQuery)$SELECTED
	
	metadata <- library1$metadata
	
	spectrum_list <- library1$sp
	
	if (nrow(metadata)==0) stop("The library does not contain the selected ID(s)!")
	
# collapse all spectrum into one data frame
	spectrum <- do.call("rbind", spectrum_list) %>% as_tibble()
	colnames(spectrum) <- c("m/z", "Intensity")
	
# prepare all metadata
	metadataExpand <- metadata[rep(1:nrow(metadata), times= unlist(lapply(spectrum_list, nrow))), ]
	
# combine spectrum and meta data	
	plotData <- bind_cols(spectrum, metadataExpand)
	
# for each ID, MSLELVE select the largest scan
	plotData %<>% 
			#.$SCANS %>% unique
			mutate(ID = readr::parse_factor(ID)) %>%
			group_by(ID, MSLEVEL) %>% filter(SCANS == max(SCANS)) %>%
			mutate(prettyIntensity = format(Intensity, digits = 4, trim = TRUE)) %>%
			#filter(Intensity > max(Intensity)*0.05) %>%
			#slice(1:100) %>%
			mutate(MSLEVEL = paste("MSLEVEL:", MSLEVEL)) %>%
			ungroup()  
	return(plotData)
}



# input_library <- library_reader("//10.157.249.83/share/mergeion/JANSSEN_COMPLETE.RData")
plotData <- dataForPlotly(input_library = input_library, 
		id = input_library$complete$metadata$ID[1:100], 
		type = "complete")
#
#library_visualizer_interactive(plotData = plotData)

#' Interactive spectral plot
#' @param plotData A data frame used to plot
#' @param x x-axis variable quote or unquoted
#' @param y y-axis variable 
#' @param coloVar color variable 
#' @param facetVar the variable used to facet the graph by row
#' @importFrom rlang ensym enexprs
#' @importFrom ggplot2 ggplot geom_segment facet_grid theme theme_bw vars
#' @importFrom plotly ggplotly
#' @author Yingjie Zhang
#' @export

library_visualizer_interactive <- function(plotData, x = "m/z", 
		y = Intensity, 
		colorVar = ID, 
		facetVar = MSLEVEL){
	
	x <- ensym(x)
	varList <- paste(enexprs(x, y, colorVar, facetVar))
	# sanity check
	condition <- sapply(   varList ,function(x){
				any(grepl(x, colnames(plotData)))
			}
	)
	
	
	
	if(any(!condition)) stop("The variables required is not found in the metadata:", names(condition)[!condition])
	
	
	p <- ggplot(data = plotData, aes(x = {{x}}, y = {{y}},
					text = paste( "Intensity:",format({{y}}, digits = 4, trim = TRUE)) ,
					color = {{colorVar}}))
	p <- p  +
			geom_segment(aes(xend =  {{x}}, yend = 0), size = 1, lineend = "butt") + 
			#geom_text_repel() +
			facet_grid(rows = vars({{facetVar}})) +
			theme_bw()  +
			theme(strip.placement = "outside",
					strip.background.x = element_rect( fill = "lightblue"),
					strip.background.y = element_rect( fill ="lightyellow"),
					panel.background = element_blank(),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					strip.text.x = element_text(size = 12),
					strip.text.y = element_text(size = 15)
			) 
	#
	ggplotly(p)
	
	
}