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
  headerPanel("SEC-SWATH-Explorer Hela CCL2 cell-cycle arrest Interphase vs. Mitosis - complex-centric viewer"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    width = 2,
    selectInput("fcolumn",
                "Choose Identifier type for gene/protein selection",
                annotations,
                selected = "Gene_names"),
    
    uiOutput("fcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
    p("Delete above entries by backspace and start typing for live search for your target protein(s)"),
    
    selectizeInput("replicate", label = "Select experimental replicate",
                   choices = c(1:3), selected = 1, multiple = FALSE),
    
    checkboxInput("split_plot", label = "Split plot by condition",
              value = TRUE),
    
    checkboxInput("logscale", "LOG10 Y-Axis", value = FALSE),
    
    checkboxInput("show_monomers", "Indicate monomer expected fractions", value = TRUE),
    
    p("Uniprot information for the current proteins is summarized in the lower table."),
    p("DISCLAIMER: THIS DATA IS UNPUBLISHED WORK - share only with permission of the creators."),
    p("Amon S., Koehler, M., Heusel M., Kutay U., Aebersold R."),
    p("Contact: heuselm@imsb.biol.ethz.ch")
    
    ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotlyOutput("plot", height = 800, width = 1600),
    dataTableOutput("table")
  )
))
