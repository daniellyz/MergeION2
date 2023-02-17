
options(shiny.maxRequestSize=10000*1024^2) 
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
library(pracma)
library(plyr)
library(tools)

shinyServer(function(input, output,clientData, session) {
  
  observeEvent(input$killButton,{
    session$reload()
  })
  
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
    params.search = list(mz_search = input$mz_search, ppm_search = input$ppm_search, rt_search = 10, rt_gap = 30)
    params.query.sp = list(prec_mz = as.numeric(input$prec_mz), use_prec = input$use_prec, polarity = input$prec_polarity, method = input$sim_methods, min_frag_match = 5, min_score = 0, reaction_type = "Metabolic")

    if (input$prec_rt!=""){
      query_expression = paste0("RT = ", input$prec_rt)
    } else {query_expression = ""}
    
    if (!is.null(query_spectrum)){

      candidates = library_query(input_library = library_file, query_expression = query_expression, 
                  query_spectrum = query_spectrum, query_library = NULL, params.search = params.search, params.query.sp = params.query.sp)

      if (is.null(candidates)){
        mms = "No candidates found!"
      } else {
        mms = "Candidates found! See next panel for annotation results!"
        candidates = candidates$consensus
    }} else {mms = "No valid query MS2 spectrum!"}
    
    return(list(candidates = candidates, mms = mms))
  })
  
  output$table1 <- renderDataTable({
    
    candidate_table=NULL
    
    withProgress({
      
      setProgress(message="Generating annotation results...")
      Sys.sleep(1)
      candidate_table = find_candidates()$candidates$metadata

      if (!is.null(candidate_table)){
            if (nrow(candidate_table)>0){
            #candidate_table = candidate_table[,c("ID", "PEPMASS", "FORMULA", "NAME","SCORE_MERGEION")]
            candidate_table = candidate_table[,c("ID", "PEPMASS", "FORMULA","SCORE_MERGEION")]
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
    
    tmp_library = library_query(library_file, query_expression = "MSLEVEL=2")

    library_visualizer(input_library = tmp_library, id = selected_id, query_spectrum = query_spectrum)
    
  })
  
  output$plot_structure <- renderPlot({
    
    selected_smi = selected_candidates()$selected_candidate$SMILES[1]
    
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,100))
    renderSMILES.rcdk(selected_smi,kekulise=TRUE)    
  })

})