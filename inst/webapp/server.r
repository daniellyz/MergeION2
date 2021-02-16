
options(shiny.maxRequestSize=100*1024^2) 
options(stringsAsFactors = F)

library(shiny)
library(V8)
library(shinyjs)
library(MergeION)
library(formattable)
library(stringr)
library(DT) 
library(prozor)
library(markdown)
library(RChemMass)

shinyServer(function(input, output,clientData, session) {
  
  observeEvent(input$killButton0,{
    session$reload()
  })
  
  observeEvent(input$killButton,{
    session$reload()
  })
  
  observeEvent(input$killButton1,{
    session$reload()
  })
  
  generate_library <-eventReactive(input$generateButton,{
    
    output_library = NULL
    mms = ""
    
    library_name = input$input_library$name
    library_file = input$input_library$datapath

    lcms_names = input$lcms_files$name
    lcms_files = input$lcms_files$datapath
    
    metadata_name = input$metadata_file$name
    metadata_file = input$metadata_file$datapath
    
    tmp_ref = NULL
    if (!is.null(metadata_file)){
      tmp_ref = read.csv(metadata_file, sep=";", dec=".",header=T)
      if (ncol(tmp_ref)==1){tmp_ref = read.csv(metadata_file, sep=",", dec=".", header=T)}  
      if (ncol(tmp_ref)==1){tmp_ref = read.csv(metadata_file, sep="\t", dec=".", header=T)}
      if ("FILENAME" %in% colnames(tmp_ref)){
        new_filenames =  lcms_files[match(tmp_ref$FILENAME, lcms_names)]
      }
      tmp_ref$FILENAME = new_filenames
      tmp_ref = tmp_ref[which(!is.na(tmp_ref$FILENAME)),,drop=FALSE]
    }

    polarity = input$polarity
    add.adduct = input$add.adduct
    mslevel = input$search_mslevel
    if (mslevel == "Both"){
      mslevel = c(1,2)
    } else {mslevel = as.numeric(mslevel)
    }
    search_mz = input$search_mz
    search_ppm = input$search_ppm
    search_rt = input$search_rt
    search_gap = input$search_gap
    
    normalized = input$process_normalized
    consensus = input$process_consensus
    baseline = input$process_baseline
    relative = input$process_relative
    max_peaks = input$process_max_peak
  
    lib_name = input$lib_name
    sample_type = input$sample_type
    user_name = input$user_name
    comments = input$user_notes
    
    if (!is.null(lcms_files) & !is.null(tmp_ref) & !is.null(lib_name)){

      output_library = try(library_generator(input_library = library_file, lcms_files = lcms_files, metadata_file = tmp_ref,
            polarity = polarity, mslevel = mslevel, add.adduct = add.adduct, processing.algorithm = "Default",
            params.search = list(mz_search = search_mz, ppm_search = search_ppm, rt_search = search_rt, rt_gap = search_gap),
            params.ms.preprocessing = list(normalized = normalized, baseline = baseline, relative = relative, max_peaks = max_peaks, recalibration = 0),
            params.consensus = list(consensus = consensus, consensus_method = "consensus", consensus_window = search_mz*2),
            params.user = list(sample_type = sample_type, user_name = user_name, comments = comments)), silent = T)

      if (class(output_library)=="try-error"){
          mms = as.character(output_library)
          output_library = NULL
      } else {mms = "Library successfully built!"}
    } else {mms = "Please check library_generator input!"} 
    
    return(list(library1 = output_library, mms = mms))
  })
  
  observeEvent(input$generateButton,{
    
    withProgress({
      setProgress(message ="Generating library...")
      Sys.sleep(1)
      mms = generate_library()$mms
    })
     output$blank_message0 <- renderPrint({mms})
  })
  
  output$downloadButton <- downloadHandler(
    
     filename = function(){
       paste0(input$lib_name, ".", input$lib_format)
     },
     content = function(con) {
       library_writer(generate_library()$library1, con)
     }
   )

  check_input <-eventReactive(input$goButton,{
    
    query_spectrum = NULL
    valid = 1
    mms = ""
    
    inFile1=input$file1
    inFile2=input$file2
    
    if (input$blank_file1==""){
        mms="Please paste your mass peaks!"
        valid = 0
     }
    
    if (is.na(as.numeric(input$prec_mz))){
        mms="Precursor mass must be a numeric value!"
        valid = 0
      } else {
        prec_mz = as.numeric(input$prec_mz)
      }
    
    if (input$blank_file1!=""){
      
      input_str = input$blank_file1
      input_str = strsplit(input_str,"\n")[[1]]
      
      input_str = lapply(input_str, function(x) strsplit(x, "\\s+")[[1]])
      input_str = lapply(input_str, function(x) x[x!="\t"])
      input_str = lapply(input_str, function(x) x[x!=""])

      if (all(sapply(input_str,length)==1)){ # One column situation
        masslist=as.numeric(unlist(input_str))
        masslist=masslist[!is.na(masslist)]
        intlist = rep(100, length(masslist))}
      
      if (any(sapply(input_str,length)>1)){ # >1 column situation
        masslist = as.numeric(sapply(input_str,function(x) x[1]))
        intlist = as.numeric(sapply(input_str,function(x) x[2]))
        valid_peaks = which(!is.na(masslist) & !is.na(intlist))
        masslist = masslist[valid_peaks, drop = FALSE]
        intlist = intlist[valid_peaks, drop = FALSE]
      }
      
      if (length(masslist)<3){
        mms = "The input spectrum must contain at least 3 fragment!"
        valid = 0
      }
    
      if (!is.null(masslist) && valid==1){ 
        query_spectrum = cbind(masslist, intlist)
      }
    }
    
    list(query_spectrum = query_spectrum, mms = mms,valid=valid)
  })

  observeEvent(input$exampleButton1, {
    fileText <- paste(readLines("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/example_cinnarizine.txt"), collapse = "\n")
    updateTextAreaInput(session, "blank_file1", value = fileText)
  })  
  
  
  observeEvent(input$goButton,{
    withProgress({
      setProgress(message = "Check data format...")
      Sys.sleep(1)
      mms = check_input()$mms
      if (check_input()$valid==1){
        setProgress(message = "Annotating by spectral library...")
        Sys.sleep(1)
        mms = find_candidates()$mms
      }
    })
    
    if (check_input()$valid==0){
      updateActionButton(session, "goButton",label = "Try again")
    } 
    
    output$blank_message1<-renderPrint({mms})
  })
  
  find_candidates <- eventReactive(input$goButton,{
  
    candidates = NULL
    mms = ""
    
    library_name = input$db_source$name
    library_file = input$db_source$datapath
    query_spectrum = check_input()$query_spectrum
    params.search = list(mz_search = input$mz_search, ppm_search = input$ppm_search, rt_search = 10, rt_gap = 0)
    params.query.sp = list(prec_mz = as.numeric(input$prec_mz), use_prec = input$use_prec, polarity = input$prec_polarity, method = input$sim_methods, min_frag_match = 6, min_score = 0, reaction_type = "Metabolic")

    if (input$prec_rt!=""){
      query_expression = paste0("RT = ", input$prec_rt)
    } else {query_expression = ""}
    
    if (!is.null(query_spectrum)){

      candidates = library_query(input_library = library_file, query_expression = query_expression, 
                  query_spectrum = query_spectrum, query_file = NULL,params.search, params.query.sp)

      if (is.null(candidates)==0){
        mms = "No candidates found!"
      } else {
        mms = "Candidates found! See next panel for annotation results!"
        candidates = candidates$consensus
    }} else {mms = "No valid query MS2 spectrum!"}
    
    return(list(candidates = candidates, mms = mms))
  })
  
  output$table0 <- renderDataTable({
    
    library_metadata=NULL

    library_name = input$db_raw$name
    library_file = input$db_raw$datapath

    if (!is.null(library_file)){
      JANSSEN_raw = library_reader(library_file)
      JANSSEN_raw = JANSSEN_raw$complete
      JANSSEN_metadata = JANSSEN_raw$metadata
      JANSSEN_sp = JANSSEN_raw$sp
      
      if (!is.null(JANSSEN_metadata)){
        if (nrow(JANSSEN_metadata)>0){
          JANSSEN_metadata = JANSSEN_metadata[,c("SAMPLEID", "ID","PEPMASS", "IONMODE", "INSTRUMENT", "SUBMITUSER")]
          library_metadata= datatable(JANSSEN_metadata, rownames = FALSE, escape= rep(TRUE, ncol(JANSSEN_metadata)), selection = "single", options = list(pageLength=20))
        }}
      return(library_metadata)
    }
  })
  
  selected_library <- eventReactive(input$table0_rows_selected,{
    
    selected_candidate = NULL    
    selected_sp = NULL    
    
    library_file = input$db_raw$datapath

    JANSSEN_raw = library_reader(library_file)
    JANSSEN_raw = JANSSEN_raw$complete
    JANSSEN_metadata = JANSSEN_raw$metadata
    JANSSEN_sp = JANSSEN_raw$sp

    selected_candidate = JANSSEN_metadata[input$table0_rows_selected,,drop=FALSE]
    selected_sp = JANSSEN_sp[[input$table0_rows_selected]]
    
    list(selected_candidate = selected_candidate, selected_sp = selected_sp)
  })
  
  output$plot_spectra <- renderPlot({
    
    selected_id = selected_library()$selected_candidate$ID[1]
    library_file = input$db_raw$datapath
    
    if (!is.null(library_file)){
      library_visualizer(input_library = library_file, id = selected_id, query_spectrum = NULL)
    }
  })
  
  output$plot_structure0 <- renderPlot({
    
    selected_smi = selected_library()$selected_candidate$SMILES[1]
    
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,100))
    renderSMILES.rcdk(selected_smi,kekulise=TRUE)    
  })

  output$table1 <- renderDataTable({
    
    candidate_table=NULL
    
    withProgress({
      
      setProgress(message="Generating annotation results...")
      Sys.sleep(1)
      candidate_table = find_candidates()$candidates$metadata

      if (!is.null(candidate_table)){
            if (nrow(candidate_table)>0){
            candidate_table = candidate_table[,c("ID", "PEPMASS", "NAME", "SCORE_MERGEION")]
            candidate_table= datatable(candidate_table, rownames = FALSE, escape= rep(TRUE, 4), selection = "single", options = list(pageLength=5))
      }}
      return(candidate_table)
    })
  })
  
  selected_candidates <- eventReactive(input$table1_rows_selected,{
    
    candidate_table = find_candidates()$candidates$metadata
    candidate_sp = find_candidates()$candidates$sp
    
    selected_candidate = candidate_table[input$table1_rows_selected,,drop=FALSE]
    selected_sp = candidate_sp[[input$table1_rows_selected]]
    
    list(selected_candidate = selected_candidate, selected_sp = selected_sp)
  })
  
  output$plot_mirror <- renderPlot({
    
    selected_id = selected_candidates()$selected_candidate$ID[1]
    library_file = input$db_source$datapath
    query_spectrum = check_input()$query_spectrum
    
    library_visualizer(input_library = library_file, id = selected_id, query_spectrum = query_spectrum)
    
  })
  
  output$plot_structure <- renderPlot({
    
    selected_smi = selected_candidates()$selected_candidate$SMILES[1]
    
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,100))
    renderSMILES.rcdk(selected_smi,kekulise=TRUE)    
  })
    
  generate_network <-eventReactive(input$networkButton,{
    
    output_network = NULL
    mms = ""
    
    library_name = input$input_library1$name
    library_file = input$input_library1$datapath
    lcms_name = input$lcms_files1$name
    lcms_file = input$lcms_files1$datapath
    polarity = input$polarity1
    
    search_mz = input$search_mz1
    search_ppm = input$search_ppm1
    search_rt = input$search_rt1
    search_gap = input$search_gap1
    
    max_peaks = input$network_max_peak
    min_frag_match = input$network_min_frag
    sim_method = input$network_similarity
    min_sim = input$network_min_score
    topK = input$network_topK
    reaction_type  = input$network_reaction_type
    
    network_name = input$network_name
    sample_type = input$network_sample_type

    if (!is.null(lcms_file) & !is.null(network_name)){
      
      output_network = try(library_generator(input_library = library_file, lcms_files = lcms_file, metadata_file = NULL,
                            polarity = polarity, mslevel = 2, add.adduct = FALSE, processing.algorithm = "Default",
                            params.search = list(mz_search = search_mz, ppm_search = search_ppm, rt_search = search_rt, rt_gap = search_gap),
                            params.ms.preprocessing = list(normalized = TRUE, baseline = 0, relative = 0.1, max_peaks = max_peaks, recalibration = 0),
                            params.consensus = list(consensus = TRUE, consensus_method = "consensus", consensus_window = search_mz*2),
                            params.network = list(network = TRUE, similarity_method = sim_method, min_frag_match = min_frag_match, min_score = min_sim, topK = topK, reaction_type = reaction_type, use_reaction = FALSE),
                            params.user = list(sample_type = sample_type, user_name = "network generator", comments = "")), silent = T)
      
      if (class(output_network)=="try-error"){
        mms = as.character(output_network)
        output_network = NULL
      } else {mms = "Network successfully built!"}
    } else {mms = "Please check library_generator input!"} 
    
    return(list(network1 = output_network, mms = mms))
  })
  
  observeEvent(input$networkButton,{
    
    withProgress({
      setProgress(message ="Generating network...")
      Sys.sleep(1)
      mms = generate_network()$mms
    })
    output$blank_message2 <- renderPrint({mms})
  })
  
  output$downloadNodes <- downloadHandler(
    
    filename = function(){
      paste0(input$network_name, "_nodes.txt")
    },
    content = function(con) {
      write.table(generate_network()$network1$network$nodes, con, quote = F, sep = "\t", col.names = T, row.names = F, dec= ".")
    }
  )

  output$downloadEdges <- downloadHandler(
    
    filename = function(){
      paste0(input$network_name, "_network.txt")
    },
    content = function(con) {
      write.table(generate_network()$network1$network$network, con, quote = F, sep = "\t", col.names = T, row.names = F, dec= ".")
    }
  )
  
  output$downloadNetworkObj <- downloadHandler(
    
    filename = function(){
      paste0(input$network_name, ".RData")
    },
    content = function(con) {
      network_obj  = generate_network()$network1
      save(network_obj, file = con)
    }
  )
})