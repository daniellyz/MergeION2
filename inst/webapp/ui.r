
options(shiny.maxRequestSize=100*1024^2) 

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

#data(JANSSEN_POS)

textInputRow<-function (inputId, label, value = ""){
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

shinyUI(navbarPage("TRACTION WEBTOOL 0.2 (Spectral library search)",
          tabPanel("A) Basic search",
                   
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
                   
          tabPanel("B) Advanced search",
                   
              shinyjs::useShinyjs(),
              shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),
                            
              column(7,  
                      br(),
                      h4("Please paste your MS/MS spectrum into the field below:"), 
                      textAreaInput("blank_file1", label = '',width=500,height=200),
                     
                      br(),
                      h4("Please choose input spectral library file: "),
                     #selectInput("db_source", label = "", choices=c("Janssen", "Drug", "GNPS", "MassBank")),  
                      fileInput("db_source", "", multiple = FALSE),
                     
                      selectInput("prec_polarity", h4("Polarity of query spectrum:"), choices=c("Positive", "Negative")),  
                     
                      br(),
                      textInput("prec_mz", h4("Precursor mass:"), value = ""),
              
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
                      selectInput("sim_methods", label = "", choices= c("Matches", "Dot", "Cosine", "Spearman", "MassBank", "NIST")),  
                    
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
          
          tabPanel("C) Annotation results",
                   tags$style("#blank_message2 {font-size:20px; color:red; display:block; }"),
       
                   tags$style("#blank_message2 {font-size:20px; color:red; display:block; }"),        
                   
                   br(),
                   div(style="display: inline-block;vertical-align:top; width: 550px;", uiOutput("blank_message2")),
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
                      plotOutput("plot_structure",width = '400px')))
))

