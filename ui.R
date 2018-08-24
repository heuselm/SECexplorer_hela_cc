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
                 p("Welcome to the Hela CCsec Viewer. Please wait a few seconds while the data is loading..."),
                 h1("Usage Instructions"),
                 p("This app has 2 main functionalities, acessible over the tabs above:"),
                 tags$ol(
                   tags$li("Viewer: In order to elucidate the interaction network complex context of proteins identified
                           in SEC, chromatogram analysis in the context of co-complex 
                           members is warranted. Therefore, SECexplorer-cc implements multi-protein selection and 
                           across-condition visualization for expert inspection"),
                   tags$li("The identified protein traces can be queried for local or global co-elution
                           in order to find potential interaction parnters.")
                 ),
                 h2("Viewer"),
                 img(src='HowtoViewer.png', align = "left", width = "100%"),
                 h3("1. Selection of ID type"),
                 tags$ul(
                   tags$li("Entry_name or Gene_Names are most informative and intuitive to search."),
                   tags$li("Searching for specific protein identifiers is then possible in field 2.")
                   ),
                 
                 h3("2. Search and selection of multiple proteins"),
                 tags$ul(
                   tags$li("Searching is achieved by deleting the current entry with backspace and starting to 
                 type. All identifiers of the type selected in (1) will be searched for the string 
                     entered, on-the-fly, with potential results showing up below the field, selectable 
                     by <Enter> or by a left mouse click."),
                   tags$li("NOTE: Individual proteins can be removed again by <Backspace>.")
                 ),
                 h3("3. Chromatogram view across conditions "),
                 tags$ul(
                   tags$li("By default, a split graph shows the protein level abundance profiles over the 
                            chromatographic fractions"),
                   tags$li("The chromatograms are displayed for one of the experimental replicates at a time"),
                   tags$li("Options for chromatogram display can be selected in (4)")
                 ),
                 h3("4. Options for chromatogram display"),
                 tags$ul(
                   tags$li("Selection of experimental replicate"),
                   tags$li("LOG10-transformation of the intensity axis to spot low-abundant protein pools"),
                   tags$li("Selection whether monomer expected fraction markers shall be displayed")
                   ),
                 
                 h3("5. Annotation table for the selected proteins "),
                 tags$ul(
                   tags$li("To help cross-referencing and establishing links which other proteins may be 
                            interesting to display in reference.")
                 ),
                 h2("Search"),
                 img(src='HowtoSearchInput.png', align = "left", width = "100%"),
                 h3("1. Experimental condition"),
                 tags$ul(tags$li("Choose the experimental condition where the search will be performed")),
                 h3("2. Base Protein"),
                 tags$ul(tags$li("Choose the protein Trace that will be searched")),
                 h3("3. Search area"),
                 tags$ul(
                   tags$li("Define in which range co-eluting proteins should be searched for"),
                   tags$li("If no selection is made all fractions will be used")
                   ),
                 h3("4. Perform the search"),
                 p(),
                 h2("Search Results"),
                 img(src='HowtoSearchResult.png', align = "left", width = "100%"),
                 h3("1. Correlation cutoff"),
                 tags$ul(
                   tags$li("Choose a cutoff for the local and global correlation of proteins to be displayed"),
                   tags$li("Global correllation: The minimum correllation across all selected SEC fractions"),
                   tags$li("Local correllation: The minimum correllation within the defined search area")
                 ),
                 h3("2. Result table"),
                 tags$ul(tags$li("A table of all proteins meeting the search criteria")),
                 h3("3. Reset"),
                 tags$ul(tags$li("Return to the search input window"))
        ),
        tabPanel('Viewer',       
                 plotlyOutput("plot", height = 600),
                 p(),
                 dataTableOutput("table")
        ),
        tabPanel('Search',
                 uiOutput("plots"),
                 fluidRow(
                   column(3, uiOutput("sliderglob")),
                   column(3, uiOutput("sliderloc"))
                 ),
                 dataTableOutput("restable"),
                 p(class = 'text-center', downloadButton('downloadData', 'Download'))
        ),
        tabPanel('DiffExpr',
                 p("Volcano plot indicating differentially behaving Proteins acrosss the cell cycle"),
                 plotlyOutput("plot_diffexpr", height = 600), 
                 verbatimTextOutput("hover")
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
