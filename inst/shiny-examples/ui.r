
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

#data(JANSSEN_POS)

textInputRow<-function (inputId, label, value = ""){
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}

shinyUI(navbarPage("meRgeION WEBTOOL 2.0 (Library Search)",
                   
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
       
          tabPanel("Submit Search",
                   
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
                      selectInput("sim_methods", label = "", choices= c("Dot", "F1", "Cosine", "HM", "MassBank", "NIST", "Entropy")),  
                    
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
                   br(),
                   plotOutput("plot_mirror",width = '1200px'),
                   br(),
                   br(),
                   plotOutput("plot_structure",width = '400px'))
))

