library(shiny)
library(ggplot2)
library(plotly)
library(data.table)

# data preparation
trace_annotation_cum<- readRDS("trace_annotation_cum.rda")

annotations <- names(trace_annotation_cum)
for (i in seq_along(annotations)){
  assign(annotations[i], trace_annotation_cum[[i]])
}

shinyUI(fluidPage(
  
  
  # Application title
  headerPanel("SEC-SWATH-Explorer Hela CCL2 cell-cycle arrest Interphase vs. Mitosis - complex-centric viewer"),
  
  # Sidebar with a slider input for number of observations
  sidebarLayout(  
    sidebarPanel(
      width = 3,
      selectInput("fcolumn",
                  "Choose Identifier type for gene/protein selection",
                  annotations,
                  selected = "Gene_names"),
      
      uiOutput("fcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
      p("Delete above entries by backspace and start typing for live search for your target protein(s)"),
      
      conditionalPanel('input.dataset === "Viewer"',
                       
                       selectizeInput("replicate", label = "Select experimental replicate",
                                      choices = c(1:3), selected = 1, multiple = FALSE),
                       checkboxInput("split_plot", label = "Split plot by condition",
                                     value = TRUE),
                       
                       checkboxInput("logscale", "LOG10 Y-Axis", value = FALSE),
                       
                       checkboxInput("show_monomers", "Indicate monomer expected fractions", value = TRUE)
                       
      ),
      conditionalPanel('input.dataset === "Search"',
                       selectizeInput("trace", label = "Select experimental condition",
                                      choices = paste0(c(rep("mit_r",3), rep("int_r", 3)),c(1,2,3)), selected = 1, multiple = FALSE),
                       uiOutput("baseProt"),
                       actionButton("search", label = "Perform Search"),
                       actionButton("reset", label = "Reset")
                       # plotOutput("plot_st_string")
      ),
      conditionalPanel('input.dataset === "String"',
                       uiOutput("stringProt"),
                       uiOutput("nrinteractors"),
                       uiOutput("confidencethr"),
                       actionButton("refresh", label = "Refresh"),
                       actionButton("paste", label = "Paste Interactors")
                       # plotOutput("plot_st_string")
      ),
      conditionalPanel('input.dataset == "Welcome"',
                       helpText("Select a tab to start")
      ),
      p("DISCLAIMER: THIS DATA IS UNPUBLISHED WORK - share only with permission of the creators."),
      p("Amon S., Koehler, M., Heusel M., Kutay U., Aebersold R."),
      p("Contact: heuselm@imsb.biol.ethz.ch")
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel('Welcome', 
                 p("Welcome to the Hela CCsec Viewer. Please wait a few seconds while the data is loading...")
        ),
        tabPanel('Viewer',       
                 plotlyOutput("plot"),
                 dataTableOutput("table")
        ),
        tabPanel('Search',
                 uiOutput("plots"),
                 fluidRow(
                   column(4, uiOutput("sliderglob")),
                   column(4, uiOutput("sliderloc")),
                   column(2, downloadButton('downloadData', 'Download'))
                 ),
                 dataTableOutput("restable")
        ),
        tabPanel('String',
                 p("String interaction partners"),
                 fluidRow(
                   column(12, align="center",
                          imageOutput("plot_string_neighbors")
                   )
                 ),
                 dataTableOutput("stringtable")
                 
        )
      )
    )
  )
))
