
options(shiny.maxRequestSize=100*1024^2) 
options(stringsAsFactors = F)

library(shiny)
library("V8")
library(shinyjs)
library(MergeION)
library(formattable)
library(stringr)
library(DT) 
library(prozor)
library(markdown)
library(RChemMass)
source("process_devco.R")

#data(JANSSEN_POS)

textInputRow<-function (inputId, label, value = ""){
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

shinyUI(navbarPage("TRACTION WEBTOOL 2.0 (Pipeline)",
                   
        tags$head(tags$style(
          HTML('
        
            #sidebar, #sidebar1 {
                  border: 1px solid black;
            }
            #sidebar, #sidebar2 {
                  border: 1px solid black;
            }
            #sidebar, #sidebar3 {
                  border: 1px solid black;
            }
            #sidebar, #sidebar4 {
                  border: 1px solid black;
            }
            #sidebar, #sidebar1bis {
                  border: 1px solid black;
            }
            #sidebar, #sidebar2bis {
                  border: 1px solid black;
            }
            #sidebar, #sidebar3bis {
                  border: 1px solid black;
            }
            #sidebar, #sidebar4bis {
                  border: 1px solid black;
            }
            body, label, input, button, select { 
            font-family: "Arial";
            }')
          )),
                   
          tabPanel("A) Library generation",
            fluidRow(
                column(3,
                sidebarLayout(
                      sidebarPanel(width = 12, id="sidebar1",
                          h3("1. Input Files"),
                          fileInput("input_library",h5("Previous spectral library file [Optional]: "), multiple = FALSE),
                          fileInput("lcms_files", h5("Converted LC-MS files (mzXML or netCDF): "), multiple = TRUE),
                          radioButtons("polarity", h5("LC-MS files polarity: "), choices = c("Positive","Negative")),
                          textInput("sample_id", h5("Devco Sample ID")),
                          fileInput("metadata_file", h5("Metadata to collect in a .csv (w/o sample id): "), multiple = FALSE),
                          checkboxInput("add.adduct", h6("Add common adducts"), FALSE, width = "200%")
                        ),
                          mainPanel(width = 0)
                        )
                  ),
                column(3,sidebarLayout(
                    sidebarPanel(width = 12, id="sidebar2",
                        h3("2. Matching Parameters"),
                        radioButtons("search_mslevel", h5("MS Level collected: "), choices = c("1","2", "Both"), selected = "2"),
                        numericInput("search_mz", h5("Mass tolerance (Da): "), min = 0, max = 0.1, value = 0.01, width = '500px'),
                        numericInput("search_ppm", h5("Mass tolerance (ppm): "), min = 0, max = 50, value = 5, width = '500px'),
                        numericInput("search_rt", h5("RT tolerance (second): "), min = 0, max = 30, value = 6, width = '500px'),
                        numericInput("search_gap", h5("RT gap (second): "), min = 0, max = 10000, value = 12, width = '500px')
                      ),
                      mainPanel(width = 0)
                ))
            ),
            
          fluidRow(
              column(3,sidebarLayout(
                sidebarPanel(width = 12, id="sidebar3",
                             h3("3. Spectrum post-processing parameters"),
                             checkboxInput("process_normalized", h5("Scaling MS peak intensity to 100"), TRUE, width = "200%"),
                             checkboxInput("process_consensus", h5("Generating consensus spectra"), TRUE, width = "200%"),
                             numericInput("process_baseline", h5("MS baseline in absolute intensity: "), min = 0, max = 10000, value = 1, width = '500px'),
                             numericInput("process_relative", h5("MS baseline in relative intensity (%): "), min = 0, max = 100, value = 0.1, width = '500px'),
                             numericInput("process_max_peak", h5("Maximum number of peaks kept: "), min = 10, max = 1000, value = 200, width = '500px')
                ),
                mainPanel(width = 0)
              )),
              column(3,sidebarLayout(
                sidebarPanel(width = 12, id="sidebar4",
                       h3("4. Library output parameters"),
                       textInput("lib_name:", h5("Name of output library:"), value = ""),
                       selectInput("lib_format", label = h5("Library file format:"), choices=c("rdata", "mgf")),
                       textInput("user_name:", h5("User initials:"), value = ""),
                       textInput("user_notes:", h5("User comments:"), value = ""),
                       selectInput("sample_type", label = h5("Sample type:"), choices=c("Impurity", "Stability", "Synthesis"))
              ),
                mainPanel(width = 0)
              )),
            column(5,
                tags$head(
                  tags$style(HTML('#generateButton{background-color:lightgreen}'))
                ),
                actionButton("generateButton", "Start Library Generation",style='padding:6px; font-size:150%'),
                br(),
                br(),
                tags$head(
                  tags$style(HTML('#downloadButton{background-color:lightblue}'))
                ),
                downloadButton("downloadButton", "Download Library",style='padding:6px; font-size:150%'),
                br(),
                br(),
                tags$head(
                  tags$style(HTML('#killButton0{background-color:orange}'))
                ),
                actionButton("killButton0", "Clear",style='padding:6px; font-size:150%'),
                br(),
                br(),
                em('Messages from the server:'),
                br(),
                textOutput("blank_message0")
              )
          )),
          
          tabPanel("B) Basic query",
                br(),
                h4("Please choose input spectral library file: "),
                fileInput("db_raw", "", multiple = FALSE),
                
                h4("Here is the list of all historical analysis"), 
                br(),
                dataTableOutput("table0"), 
                br(),
                column(8,  
                       br(),
                       plotOutput("plot_spectra",width = '1200px')),
                column(3,
                       br(),
                       plotOutput("plot_structure0",width = '400px'))
          
          ),
                   
          tabPanel("C) Library search",
                   
              #shinyjs::useShinyjs(),
              #shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),
                            
              column(7,  
                      br(),
                      h4("Please paste your MS/MS spectrum into the field below:"), 
                      textAreaInput("blank_file1", label = '',width=500,height=200),
                     
                      br(),
                      h4("Please choose input spectral library file: "),
                      fileInput("db_source", "", multiple = FALSE),
                     
                      selectInput("prec_polarity", h4("Polarity of query spectrum:"), choices=c("Positive", "Negative")),  
                     
                      br(),
                      textInput("prec_mz", h4("Precursor mass:"), value = "10000"),
              
                      br(),
                      textInput("prec_rt", h4("Retention time in min [Optional]:"), value = "")),

                            
               column(4,
                      br(),
                      checkboxInput("use_prec", h4("Searching precursor mass"), TRUE, width = "200%"),
                      
                      br(),
                      numericInput("ppm_search", h5("Mass tolerance for precursor mass (ppm)"), min = 0, max = 20, value = 3, width = '500px'),
                      
                      br(),
                      h4("Mass tolerance for fragment matching (Da):"),
                      numericInput("mz_search", "", min = 0, max = 0.05, value = 0.01, width = '500px'),
                     
                      br(),
                      h4("Spectral similarity metrics:"),
                      selectInput("sim_methods", label = "", choices= c("Precision", "Recall", "F1", "Cosine", "Spearman", "MassBank", "NIST")),  
                    
                      br(),
                      tags$head(
                        tags$style(HTML('#exampleButton1{background-color:lightblue}'))
                      ),
                      actionButton("exampleButton1", "Load example: Cinnarizine",style='padding:6px; font-size:120%'),
                      br(),
                      
                      br(),
                      tags$head(
                        tags$style(HTML('#goButton{background-color:lightgreen}'))
                      ),
                      actionButton("goButton", "Submit",style='padding:6px; font-size:150%'),
                      br(),
                      
                      br(),
                      tags$head(
                        tags$style(HTML('#killButton{background-color:orange}'))
                      ),
                      actionButton("killButton", "Clear",style='padding:6px; font-size:150%'),
                      br(),  
                      
                      br(),
                      em('Messages from the server:'),
                      br(),
                      textOutput("blank_message1")
                  )),
          
          tabPanel("Annotations",
          
                   br(),
                   br(),     
                   h4("Here is the list of annotated candidates"), 
                   br(),
                   dataTableOutput("table1"), 
                   br(),
                   column(8,  
                      br(),
                      plotOutput("plot_mirror",width = '1200px')),
                   column(3,
                      br(),
                      plotOutput("plot_structure",width = '400px'))),
        
        
        tabPanel("D) Molecular networking",
                 fluidRow(
                   column(3,
                          sidebarLayout(
                            sidebarPanel(width = 12, id="sidebar1bis",
                                         h3("1. Input Files"),
                                         fileInput("input_library1",h5("Spectral library file [Optional]: "), multiple = FALSE, width = '500px'),
                                         fileInput("lcms_files1", h5("Converted LC-MS files (mzXML or netCDF): "), multiple = TRUE, width = '500px'),
                                         radioButtons("polarity1", h5("LC-MS files polarity: "), choices = c("Positive","Negative")),
                                         fileInput("metadata_file1", h5("Features to collect in a .csv metadata file: "), multiple = FALSE)
                            ),
                            mainPanel(width = 0)
                          )
                   ),
                   column(3,sidebarLayout(
                     sidebarPanel(width = 12, id="sidebar2bis",
                                  h3("2. Searching Parameters"),
                                  numericInput("search_mz1", h5("Mass tolerance (Da): "), min = 0, max = 0.1, value = 0.01, width = '500px'),
                                  numericInput("search_ppm1", h5("Mass tolerance (ppm): "), min = 0, max = 50, value = 5, width = '500px'),
                                  numericInput("search_rt1", h5("RT tolerance (second): "), min = 0, max = 30, value = 6, width = '500px'),
                                  numericInput("search_gap1", h5("RT gap (second): "), min = 0, max = 10000, value = 12, width = '500px'),
                                  numericInput("search_baseline1", h5("Baseline: "), min = 0, max = 100000, value = 25000, width = '500px')
                    ),
                     mainPanel(width = 0)
                   ))),
                   
                fluidRow(
                   column(3,sidebarLayout(
                     sidebarPanel(width = 12, id="sidebar3bis",
                                 h3("3. Networking parameters"),
                                 numericInput("network_max_peak", h5("Maximum number of peaks kept: "), min = 10, max = 1000, value = 200, width = '500px'),
                                 numericInput("network_min_frag", h5("Minimum of fragment matches: "), min = 5, max = 1000, value = 8, width = '500px'),
                                 selectInput("network_similarity", h5("Spectral similarity type: "), choices= c("Messar", "Precision", "Recall", "F1", "Cosine", "Spearman", "MassBank", "NIST"), width = '500px'),
                                 numericInput("network_min_score", h5("Minimum similarity to connect two nodes: "), min = 0, max = 1, value = 0.6, width = '500px'),
                                 numericInput("network_topK", h5("Each other's top K most similar nodes: "), min = 5, max = 50, value = 10, width = '500px'),                                 
                                 radioButtons("network_reaction_type", h5("Precursor mass difference type: "), choices = c("Chemical","Metabolic"))
                     ),
                     mainPanel(width = 0)
                   )),

                  column(3,sidebarLayout(
                   sidebarPanel(width = 10, id="sidebar4bis",
                                h3("4. Network output parameters"),
                                textInput("network_name:", h5("Name of output library:"), value = ""),
                                selectInput("network_sample_type", label = h5("Sample type:"), choices=c("Impurity", "Stability", "Synthesis", "Metabolomics"))
                   ),
                   mainPanel(width = 0)
                 )),
                
                  column(5,
                        tags$head(
                          tags$style(HTML('#networkButton{background-color:lightgreen}'))
                        ),
                        actionButton("networkButton", "Start Network Generation",style='padding:6px; font-size:150%'),
                        br(),
                        br(),
                        tags$head(
                          tags$style(HTML('#downloadNodes{background-color:lightblue}'))
                        ),
                        downloadButton("downloadNodes", "Download Nodes as txt",style='padding:6px; font-size:150%'),
                        br(),
                        br(),
                        tags$head(
                          tags$style(HTML('#downloadEdges{background-color:lightblue}'))
                        ),
                        downloadButton("downloadEdges", "Download Edges as txt",style='padding:6px; font-size:150%'),
                        br(),
                        br(),
                        tags$head(
                          tags$style(HTML('#downloadNetworkObj{background-color:lightblue}'))
                        ),
                        downloadButton("downloadNetworkObj", "Download Network Object",style='padding:6px; font-size:150%'),
                        br(),
                        br(),
                        tags$head(
                          tags$style(HTML('#killButton1{background-color:orange}'))
                        ),
                        actionButton("killButton1", "Clear",style='padding:6px; font-size:150%'),
                        br(),
                        br(),
                        em('Messages from the server:'),
                        br(),
                        textOutput("blank_message2")
                   )
                 ))
))

