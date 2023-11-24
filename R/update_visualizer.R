

#' Function to generate data for interactive visulization
#' @inheritParams library_visualizer
#' @param t numeric value specifying the tolerance used to align the m/z values of the two spectra (used in calculating consine similarity score), 0.25 by default. 
#' @param b numeric value specifying the baseline threshold for peak identification. Expressed as a percent of the maximum intensity (used in calculating consine similarity score), 10 by default.
#' @importFrom rex rex
#' @importFrom dplyr filter group_by bind_cols bind_rows mutate ungroup
#' @importFrom magrittr %>%  %<>%
#' @importFrom tidyr pivot_wider
#' @importFrom tibble as_tibble
#' @return a dataframe with m/z, intensity, ID and all the metadata columns from \code{input_library} metadata
#' @author Yingjie Zhang
#' @export 

#load(file.path("//10.157.249.83/share/mergeion","JANSSEN_COMPLETE.RData")) 
#input_library = JANSSEN_new
#query_spectrum <- JANSSEN_new$complete$sp[[1]][1:100,]
#id = input_library$complete$metadata$ID[1]
# type = "complete"
#testData <- dataForPlotly(input_library = JANSSEN_new, id = JANSSEN_new$complete$metadata$ID[1], query_spectrum =query_spectrum )

dataForPlotly <- function(input_library, 
		id = NULL,
		type = c("complete", "consensus"),
		query_spectrum = NULL,
		t = 0.25, b = 10){
	# input_library is already filtered, plot will be based on the entir input_library
	
	stopifnot(!missing(input_library))	
	
	input_library <-  library_reader(input_library)
	
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
	if(!is.null(id)){
		idsQuery <- paste0("ID@xxx@", 
				paste0(sapply(id, function(x){paste0("(^",rex(x),"$)")}), collapse = "|")
		)  
		
		library1 <- process_query(library0 = input_librarySub, query = idsQuery)$SELECTED
		
	}else{
		idsQuery = ""
		library1 <- input_librarySub
	}
	
	
	
	metadata <- library1$metadata
	
	spectrum_list <- library1$sp
	
	if (nrow(metadata)==0) stop("The library does not contain the selected ID(s)!")
	
# collapse all spectrum into one data frame
	spectrum <- do.call("rbind", spectrum_list) %>% as_tibble()
	
	if(ncol(spectrum) == 2) colnames(spectrum) <- c("m/z", "Intensity")
	
# prepare all metadata
	metadataExpand <- metadata[rep(1:nrow(metadata), times= unlist(lapply(spectrum_list, nrow))), ]
	
# combine spectrum and meta data	
	if(is.null(	query_spectrum)){
		plotData <- bind_cols(spectrum, metadataExpand)
		
 if(!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
# for each ID, MSLELVE select the largest scan
		
if(	type == "complete" & all(is.na(plotData$SCANS))) stop("No SCANS found in the data, we don't know which one to select")

		plotData %<>% 
				#.$SCANS %>% unique
				mutate(ID = readr::parse_factor(ID)) %>%
		    filter(SCANS == max(SCANS)) %>%
				group_by(ID, MSLEVEL) 

		
		plotData %<>%
				#mutate(prettyIntensity = format(Intensity, digits = 4, trim = TRUE)) %>%
				#filter(Intensity > max(Intensity)*0.05) %>%
				#slice(1:100) %>%
				mutate(MSLEVEL = paste("MSLEVEL:", MSLEVEL)) %>%
				ungroup() }else{
		
		# generate mirror plot data
		
		if (ncol(query_spectrum)<2){
			stop("Spectrum must have 2 columns m/z and intensity!")
		}
		
		colnames(query_spectrum) <- colnames(spectrum) 
#		plotData <- bind_rows(spectrum, as_tibble(query_spectrum), .id = "set") %>%
#				 mutate(set = recode(set,  `1` = "top", `2`= "bottom"))
		
		## this is based on SpectrumSimilarity
		spec.top <- spectrum
		spec.bottom <- as_tibble(query_spectrum)
		top_tmp <- data.frame(mz = spec.top$`m/z`, intensity = spec.top$Intensity)
		top_tmp$normalized <- round((top_tmp$intensity/max(top_tmp$intensity)) * 100)
		#top_tmp <- subset(top_tmp, top_tmp$mz >= xlim[1] & top_tmp$mz <= xlim[2])
		top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)
		top <- subset(top_plot, top_plot$intensity >= b)
		bottom_tmp <- data.frame(mz = spec.bottom$`m/z`, intensity = spec.bottom$Intensity)
		bottom_tmp$normalized <- round((bottom_tmp$intensity/max(bottom_tmp$intensity)) * 100)
		#bottom_tmp <- subset(bottom_tmp, bottom_tmp$mz >= xlim[1] & bottom_tmp$mz <= xlim[2])
		bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)
		bottom <- subset(bottom_plot, bottom_plot$intensity >= b)
		
		# perform alignment
		for (i in 1:nrow(bottom)) top[, 1][bottom[, 1][i] >= top[, 1] - t & bottom[, 1][i] <= top[, 1] + 
							t] <- bottom[, 1][i]
		alignment <- merge(top, bottom, by = 1, all = TRUE)
		if (length(unique(alignment[, 1])) != length(alignment[, 1])) 
			warning("the m/z tolerance is set too high")
		alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0
		names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
		
#	if (x.threshold < 0) 
#		stop("x.threshold argument must be zero or a positive number")
#	alignment <- alignment[alignment[, 1] >= x.threshold, ]
		u <- alignment[, 2]
		v <- alignment[, 3]
		similarity_score <- as.vector((u %*% v)/(sqrt(sum(u^2)) * sqrt(sum(v^2))))
		
		# plot using normalized but unfiltered data
		plotData <- bind_rows(top_plot, bottom_plot, .id = "set") %>% mutate(set = recode(set,  `1` = "Reference", `2`= "Query spectrum"))
		
		plotData$similarity <- 	similarity_score
		plotData <- plotData[,c(2,3,4,1)]
		plotData$ID <- id
		colnames(plotData) <-  c( "m/z", "Intensity", "similarity","set", "ID")
		
	} 
	
	return(plotData)
}



#' Interactive spectral plot
#' @param plotData A data frame used to plot
#' @param x x-axis variable quote or unquoted or column index
#' @param y y-axis variable or column index
#' @param coloVar color variable or column index
#' @param facetVar the variable used to facet the graph by row or column index
#' @importFrom rlang ensym enexprs
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom rlang as_string ensym enexpr
#' @importFrom magrittr %<>%
#' @author Yingjie Zhang
#' @export


# testPlot <- library_visualizer_interactive(plotData = testData)

library_visualizer_interactive <- function(plotData, x = "m/z", 
		y = "Intensity", 
		colorVar = "ID", 
		facetVar = "MSLEVEL",
		title = "",
		returnData = FALSE){
	
	
	
	
	if(all(colnames(plotData) ==  c("m/z", "Intensity", "similarity","set", "ID"))){
		# here begins mirror plot, none of the parameter are needed
		plotData %<>% mutate(across(.cols = -c(set,ID), ~round(.x, digits = 4)))
#		plotDataTop <- plotData %>% filter(set == "top", Intensity != 0 ) %>% dplyr::select(-set)
#		plotDataBottom <- plotData %>% filter(set == "bottom", Intensity != 0) %>% dplyr::select(-set)
#	
		plotData$key <- 1:nrow(plotData)
		
		plotData %<>% filter(Intensity != 0) %>% 
				mutate(Intensity2 = ifelse(set == "Reference", -Intensity, Intensity))
		
		p <- ggplot(data = plotData,aes(x = `m/z`, y = 0, text = paste0("Intensity: ",Intensity, "%\n",
										"m/z:", `m/z`, "\n", set), key = key)) + 
				geom_segment(aes( xend =  `m/z`, yend = Intensity2, color = set), size = 1, lineend = "butt", show.legend = FALSE) +
				theme_minimal() +
				scale_color_manual(values = c("Reference" = "darkgray", "Query spectrum" = "green4"))+
				scale_y_continuous(labels = abs) +
				ylab("Intensity (%)") + 
				annotate("text", x = quantile(plotData$`m/z`, 0.4), y = 105, size = 5,
						label = paste("Cosine similarity:", unique(plotData$similarity)), color = "antiquewhite4")
		
	}else{
		# if the variable names are column index
		if(is.numeric(enexpr(x))) x <- colnames(plotData)[x]
		if(is.numeric(enexpr(y))) y <- colnames(plotData)[y]
		if(is.numeric(enexpr(colorVar))) colorVar <- colnames(plotData)[colorVar]
		if(is.numeric(enexpr(facetVar))) facetVar <- colnames(plotData)[facetVar]
		
		varList <- paste(enexprs(x, y, colorVar, facetVar))
		# sanity check
		condition <- sapply(   varList ,function(x){
					any(grepl(x, colnames(plotData)))
				}
		)
		
		
		
		if(any(!condition)) stop("The variables required is not found in the metadata:", names(condition)[!condition])
		
		
		#############
		# turn varaible name to symbol
		
		x <- ensym(x)
		y <- ensym(y)
		colorVar <- ensym(colorVar)
		facetVar <- ensym(facetVar)
		
		# chcek if Frag_Sub, Frag_Formula, Nloss, Nloss_Formula present
		plotData$key <- 1:nrow(plotData)
		if(all(c("Frag_Sub", "Frag_Formula", "Nloss", "Nloss_Formula") %in% colnames(plotData))){
			p <- ggplot(data = plotData, aes(x = {{x}}, y = {{y}}, key = key,
							text = paste( as_string(y), ": ",format({{y}}, digits = 4, trim = TRUE),
									"\n",as_string(x), ": ", format({{x}}, digits = 4, trim = TRUE),
									"\n ",as_string(colorVar), ": ", {{colorVar}},
									"\n Frag_Sub: ",   Frag_Sub, 
									"\n Frag_Formula: ", Frag_Formula, 
									"\n Nloss: ", format(Nloss, digits = 4, trim = TRUE),
									"\n Nloss_Formula: ", format(Nloss_Formula, digits = 4, trim = TRUE),
									sep = "") ,
							color = {{colorVar}}))
			
		}else{
			
			
			
			p <- ggplot(data = plotData, aes(x = {{x}}, y = {{y}},
							text = paste( as_string(y), ": ",format({{y}}, digits = 4, trim = TRUE),
									"\n",as_string(x), ": ", format({{x}}, digits = 4, trim = TRUE),
									"\n ",as_string(colorVar), ": ", {{colorVar}},
									sep = "") ,
							color = {{colorVar}}))
		}
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
				) }
	
	
	p <- p + ggtitle(title)
	if(returnData) return(list(data = plotData, plot = p)) else return(p)
	

}


