####################################################################################
# SECexplorer-cc cell cycle complex association dynamics viewer ####################
# Part 1: UI.R #####################################################################
####################################################################################
# Domain: https://sec-explorer.shinyapps.io/hela_cellcycle/
# Authors: Max Frank, Moritz Heusel
# ##############################################################
# About:
# Browsing dynamic complex association maps in SECexplorer-cc:
# SEC-SWATH-MS portrays the process of mitosis from the angle of protein
# mass re-distribution across differently sized stable complexes resolved by SEC.
# There remains a lot to be learned from the rich dataset generated, with the
# prospect of discovering new proteins involved in cell cycle regulation and
# of identifying new components of both static as well as dynamic assemblies.
# To support community-based mining and interpretation of our dataset, we provide
# a web tool, SECexplorer-cc (https://sec-explorer.shinyapps.io/hela_cellcycle/)
# which offers four functionalities:
# (i) Interactive viewing of protein SEC fractionation profiles in interphase and mitosis
# (ii) Search for locally co-eluting proteins to identify putative new binding partners
# showing strong co-elution within a certain range of target protein elution
# (iii) Interactive display and protein selection from the differential association
# score map. 
# (iv) Display of one or multiple protein's fractionation profiles in reference to
# the profiles of immediate interaction and/or functional partners dynamically
# retrieved from StringDB (Szklarczyk et al., 2017). We expect that SECexplorer-cc
# will support a community effort to fully leverage the rich information encoded by
# the mitotic proteome rearrangement SEC-SWATH-MS data, which, ideally,
# will support better understanding of modular proteome function along cell division.
#####################################################################################

## prepare environment
if (!require("shiny")){
  install.packages("shiny")
}
if (!require("shinythemes")){
  install.packages("shinythemes")
}
if (!require("ggplot2")){
  install.packages("ggplot2")
}
if (!require("plotly")){
  install.packages("plotly")
}
if (!require("data.table")){
  install.packages("data.table")
}
if (!require("DT")){
  install.packages("DT")
}
if (!require("shinythemes")){
  install.packages("shinythemes")
}

# load packages
library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(data.table)
library(DT)


## prepare data
idcols <- readRDS("www/data/idcols.rda")

## define user interface
########################

shinyUI(fluidPage(
  
  # Theme
  theme = shinytheme("cosmo"),

  # Application title
  headerPanel("SECexplorer-cc: Browsing dynamic HeLa CCL2 complex association states in interphase vs. mitosis", windowTitle = "SECexplorer-cc"),

  # Sidebar with a slider input for number of observations
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # passwordInput("pwd", "Enter Password", value = "", width = NULL),
      # actionButton("enterpwd", "Enter"),
      # verbatimTextOutput("pwdfeedback"),
      p(),
      selectInput(inputId = "fcolumn",
                  label = "Choose protein profiles based on:",
                  choices = idcols,
                  selected = "Gene_names"),

      uiOutput("fcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
      p("Delete above entries by backspace and start typing for live search for your target protein(s)"),

      conditionalPanel('input.dataset === "(i) Protein profile viewer"',

                       ## selectizeInput("replicate", label = "Select experimental replicate",
                       ##                choices = c(1:3), selected = 1, multiple = FALSE),
                       checkboxInput("split_plot", label = "Split plot by condition",
                                     value = TRUE),
                       checkboxInput("logscale", "LOG10 Y-Axis", value = FALSE),
                       checkboxInput("show_monomers", "Indicate monomer expected fractions", value = TRUE),
                       checkboxInput("error_bars", "Plot error area", value = TRUE),
                       uiOutput("errortype")
      ),
      conditionalPanel('input.dataset === "(ii) Search for co-eluting proteins"',
                       selectizeInput("trace", label = "Select experimental condition",
                                      choices= c("Mitosis", "Interphase")),
                                      ## choices = paste0(c(rep("mit_r",3), rep("int_r", 3)),c(1,2,3)), selected = 1, multiple = FALSE),
                       uiOutput("baseProt"),
                       actionButton("search", label = "Perform Search"),
                       actionButton("pastediff", label = "Paste Selection"),
                       actionButton("reset", label = "Reset")
                       # plotOutput("plot_st_string")
      ),
      conditionalPanel('input.dataset === "(iii) Browse differential association map"',
                       p("To select Proteins of interest click the protein or drag a selection box around multiple proteins")
                       # plotOutput("plot_st_string")
      ),
      conditionalPanel('input.dataset === "(iv) Query StringDB partners"',
                       uiOutput("stringProt"),
                       uiOutput("nrinteractors"),
                       uiOutput("confidencethr"),
                       actionButton("refresh", label = "Refresh"),
                       actionButton("paste", label = "Paste Interactors")
                       # plotOutput("plot_st_string")
      ),
      br(),
      br(),
      br(),
      p("Original study: A global screen for assembly state changes of the mitotic proteome by SEC-SWATH-MS, Cell Systems, 02/2020"),
      p("Authors:"),
      p("Heusel M, Frank M, Koehler M, Amon S, Frommelt F, Rosenberger G, Bludau I, Aulakh S, Linder MI, Liu Y, Collins BC, Gstaiger M, Kutay U, Aebersold R"),
      p("Contact: moritz.heusel@med.lu.se, aebersold@imsb.biol.ethz.ch")
      
    ),

    # The panel for the plot output
    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel('Welcome & Instructions',
                 h1("Welcome to SECexplorer-cc"),
                 p("SECexplorer-cc offers a dynamic interface for the investigator-driven exploration of a mass spectrometric dataset describing assembly state changes in the mitotic proteome of the HeLa CCL2 cell line, chemically arrested in interphase and prometaphase/Mitosis, followed by fractionation of intact complexes and quantification of component proteins by SEC-SWATH-MS (Heusel et al., 2020). The dataset contains complex patterns of quantitative and qualitative changes in protein assemblies that warrant further, investigator-driven exploration, including the discovery of additional proteins involved in cell cycle regulation or the identification of new components of invariant or state specific assemblies. The SECexplorer-cc interface supports manual review and community-based mining and interpretation of our dataset by four core functionalities: "),
                 tags$ol(
                   tags$ul("(i) Protein profile viewer: Allows interactive viewing of protein SEC fractionation profiles in interphase and mitosis and download of chromatogram graphs to investigate the elution behavior of customized protein sets of interest."),
                   tags$ul("(ii) Search for co-eluting proteins: Identify putative binding partners of proteins based on local and global elution profile similarity."),
                   tags$ul("(iii) Browse differential association map: View differential scores for selected proteins or select proteins based on their scores in difference testing. Are there interesting outliers for which no profile change would have been anticipated?"),
                   tags$ul("(iv) Query StringDB partners: Retrieve interaction or functional partner proteins' elution profiles for contextual interpretation. Do the data support hypotheses on re-arrangements among custom complexes not investigated in the original study?")
                 ),
                 h1("Usage instructions:"),
                 h2("(i) Protein profile viewer"),
                 img(src='1_Viewer.png', align = "left", width = "100%"),
                 h3("1. Selection of ID type"),
                 tags$ul(
                   tags$li("Entry_name or Gene_Names are most informative and intuitive to search."),
                   tags$li("Searching for specific gene/protein identifiers is then possible in field 2.")
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
                            chromatographic fractions in either condition, interphase and mitosis."),
                   tags$li("The chromatograms are displayed as mean of the three replicates, with shaded areas indicating the standard error of the mean (SEM)"),
                   tags$li("Options for chromatogram display can be selected in (4)")
                 ),
                 h3("4. Options for chromatogram display"),
                 tags$ul(
                   tags$li("Selection whether the plot should be split by condition, or if conditions shall be encoded as line type."),
                   tags$li("Option for LOG10-transformation of the intensity axis to spot low-abundant protein pools"),
                   tags$li("Selection whether monomer expected fraction markers shall be displayed"),
                   tags$li("Option whether (and which type of) error areas shall be indicated. Works best on split plot.")
                   ),
                 h3("Note that the plots are interactive (thanks, plotly!)"),
                 p(),
                 
                 # h3("5. Annotation table for the selected proteins "),
                 # tags$ul(
                 #   tags$li("To help cross-referencing and establishing links which other proteins may be 
                 #            interesting to display in reference.")
                 # ),
                 h2("(ii) Search for co-eluting proteins"),
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
                 h2("(ii) Search for co-eluting proteins: Search Results"),
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
                 tags$ul(tags$li("Return to the search input window")),
                 p(),
                 h2("(iii) Browse differential association map"),
                 img(src='3_BrowseScoreMap.png', align = "left", width = "100%"),
                 h3("1. Select proteins based on scores"),
                 tags$ul(
                   tags$li("Each point represents one protein feature's differential scores between the conditions. Click the point or draw a box around the scores to select a set of proteins."),
                   tags$li("Proteins are added to those already selected via the viewer. Visit this page to see the differential scores of a selected protein set of interest.")),
                 h3("2. Conditional SEC profiles of the selected proteins. For more visualization options, visit the viewer tab (i)"),
                 p(),
                 h2("(iv) Query StringDB partners"),
                 img(src='4_QueryStringPartners.png', align = "left", width = "100%"),
                 h3("1. Selected proteins"),
                 h3("2. Select one protein as basis for the query to StringDB"),
                 h3("3. Select parameters for the query to StringDB"),
                 h3("4. Add the retrieved partner proteins to the current analysis"),
                 h3("5. StringDB network view of the current query"),
                 h3("6. Preview of the SEC-SWATH-MS chromatograms of the proteins of the current query (MS-detectable subset)")
                 
        ),
        tabPanel('(i) Protein profile viewer',       
                 plotlyOutput("plot", height = 600),
                 p(),
                 DT::dataTableOutput("table"),
                 downloadButton("downloadPlot", "Download as PDF")
        ),
        tabPanel('(ii) Search for co-eluting proteins',
                 uiOutput("plots"),
                 fluidRow(
                   column(3, uiOutput("sliderglob")),
                   column(3, uiOutput("sliderloc"))
                 ),
                 DT::dataTableOutput("restable"),
                 p(class = 'text-center', uiOutput("downloadSemiTargetedRes"))
        ),
        tabPanel('(iii) Browse differential association map',
                 ## p("Volcano plot indicating differentially behaving Proteins acrosss the cell cycle"),
                 plotlyOutput("plot_diffexpr", height = 400),
                 ## verbatimTextOutput("click"),
                 plotlyOutput("volcanotraces", height = 400)
                 # verbatimTextOutput("hover"),
                 # verbatimTextOutput("diffexprsel")
        ),
        tabPanel('(iv) Query StringDB partners',
                 p("String interaction partners"),
                 fluidRow(
                   column(12, align="center",
                          imageOutput("plot_string_neighbors")
                   )
                 ),
                 plotlyOutput("stringtraces", height=400),
                 DT::dataTableOutput("stringtable")
        )
      )
    )
  )
))
