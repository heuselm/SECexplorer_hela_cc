library(shiny)
library(ggplot2)
library(plotly)
library(data.table)

# data preparation
load("data.rda")

annotations <- names(trace_annotation_cum)
for (i in seq_along(annotations)){
  assign(annotations[i], trace_annotation_cum[[i]])
}

shinyUI(pageWithSidebar(

  
  # Application title
  headerPanel("SEC-SWATH-Explorer Hela CCL2 cell-cycle arrest experiment Interphase vs. Mitosis - complex-centric"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    width = 2,
    selectInput("fcolumn",
                "Choose Annotation Column for filtering",
                annotations,
                selected = "Gene_names"),
    
    uiOutput("fcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
    p("Delete above entry and start typing for live search for your target protein"),
    
    selectizeInput("replicate", label = "Select replicate for peptide level plot",
                   choices = c(1:3), selected = 1, multiple = FALSE),
    
    p("Select annotation column and value, and replicate to display matching SEC-SWATH chromatograms"),
    
    checkboxInput("split_plot", label = "Split plot by condition",
              value = TRUE),
    
    checkboxInput("logscale", "LOG10 Y-Axis", value = FALSE),
    
    checkboxInput("show_monomers", "Indicate monomer expected fraction", value = TRUE),
    
    p("Uniprot information for the current proteins is summarized in the lower table."),
    p("DISCLAIMER: THIS DATA IS UNPUBLISHED WORK - share only with permission of the creators."),
    p("Amon S., Koehler, M., Heusel M., Kutay U., Aebersold R."),
    p("Contact: heuselm@imsb.biol.ethz.ch")
    
    ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotlyOutput("plot"),
    dataTableOutput("table")
  )
))
