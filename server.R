####################################################################################
# SECexplorer-cc cell cycle complex association dynamics viewer ####################
# Part 2: Server.R #################################################################
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
if (!require("STRINGdb")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("STRINGdb")
}
if (!require("shiny")){
  install.packages("shiny")
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

# load packages
library(shiny)
library(ggplot2)
library(plotly)
library(data.table)
library(DT)
library(STRINGdb)

# Source custom functions
source("searchSemiTargeted.R")
source("tracesMethods.R")
source("stringMethods.R")


## prepare data
setwd("www/data")
# pass <- readRDS("pass.rda")
## trace_annotation_cum <- readRDS("trace_ann.rda")
calibration_functions <- readRDS("calibration_functions.rda")
up <- readRDS("uniprotMapping.rda")
## load("data_.rda")
tr <- readRDS("proteinTracesLong_mean_sd_sem.rda")
stringLinks <- fread("9606.protein.links.v10.5.HeLaSubset.txt")
stringIdMap <- readRDS("stringIdMapUniq.rda")
#Load the differential expression data
# diffExprProt <- readRDS("differentiallyExpressedProteins.rda")
# diffExprProt <- merge(diffExprProt, trace_annotation_cum, by.x = "feature_id", by.y = "protein_id")
# saveRDS(diffExprProt, "differentiallyExpressedProteinsAnn.rda")
diffExprProt <- readRDS("differentiallyExpressedProteinsAnn.rda")
setwd("../../")
default_proteins <- c("NUP54", "NUP62", "NUP58 KIAA0410 NUPL1")

## define server roles
#######################

shinyServer(function(input, output, session) {

  ## Set some global variables
  apex <- reactiveValues(bound_left = NULL, bound_right = NULL)
  searchResFilt <- reactiveVal(NULL)

  ############################
  ## Viewer Tab              #
  ############################

  ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
  output$fcolumnvalues <- renderUI({
    values <- sort(unique(up[[input$fcolumn]]))
    # values <- values[nchar(values)>0]
    selectizeInput("fvalue", "Search and select proteins of interest", values,
                   multiple = TRUE, options = list(maxOptions = 6000),
                   selected = default_proteins)
  })
  output$errortype <- renderUI({
    if(input$error_bars){
      selectizeInput("ertype",
                     label=NULL,
                     choices=c("Standard error of the mean (SEM)", "Standard deviation (SD)"),
                     selected=2)
    }else{
      NULL
    }
  })

  ## generate selected protein SEC traces plot
  lx.frc <- seq(5,(max(tr$fraction)-1),5)
  lx <- paste(lx.frc , round(calibration_functions$FractionToMW(lx.frc), 1) , sep = '(' )
  lx <- paste(lx, "", sep = ')' )

  # collect chromatograms for selection
  target_id_traces <- eventReactive(input$fvalue,{
    # apex$bound_right <- NULL
    # apex$bound_left <- NULL

    target_id <- up[which(up[[input$fcolumn]] %in% input$fvalue),
                         unique(Entry)]
    target_id_traces <- tr[id %in% target_id]

    target_id_traces[, monomer_fraction:=calibration_functions$MWtoFraction(protein_mw)]
    target_id_traces
  })

  ## Plot the selected traces
  viewerPlot <- function(traces, showmonomers=T, log=F, split=T, errbars=T, errtype="Standard deviation (SD)"){
    # PLOT
    p <- ggplot(traces, aes(x=fraction, y=intensity_mean + 1)) +
      xlab("SEC fraction number(apparent MW[kDa])") +
      ylab("Protein level SWATH-MS intensity (top2 peptide sum)")

    if (split){
      p <- p + geom_line(aes(group = Gene_names, color = Gene_names), size = 1, alpha = 0.8) +
        theme_bw() +
        ## ggtitle(unique(target_id_traces()$Gene_names)) +
        scale_x_continuous(breaks = lx.frc, labels = lx)  +
        facet_wrap(~condition, ncol = 1)
    } else{
      p <- p + geom_line(aes(color = Gene_names, linetype = condition), size = 1, alpha = 0.8) +
        theme_bw() +
        ## ggtitle(unique(target_id_traces()$Gene_names)) +
        scale_x_continuous(breaks = lx.frc, labels = lx)
    }
    if(errbars){
      if(input$ertype == "Standard deviation (SD)"){
        p <- p + p +
          geom_ribbon(aes(ymin=intensity_mean-intensity_sd,
                          ymax = intensity_mean+intensity_sd,
                          fill = Gene_names),
                      alpha = 0.2)
      }
      if(input$ertype == "Standard error of the mean (SEM)"){
        p <- p +
          geom_ribbon(aes(ymin=intensity_mean-intensity_se,
                            ymax = intensity_mean+intensity_se,
                            fill = Gene_names),
                        alpha = 0.2)
      }
    }
    if(showmonomers){
      p <- p + geom_point(aes(x = monomer_fraction, color = Gene_names, y = 0), shape = 23, fill = "white", size = 3)
    }

    if (log){
      p <- p + scale_y_log10()
    }
    return(p)
  }

  output$plot <- renderPlotly({
    vplot <<- viewerPlot(target_id_traces(),
                         log=input$logscale,
                         showmonomers=input$show_monomers,
                         split=input$split_plot,
                         errbars=input$error_bars,
                         errtype=input$ertype)
    ggplotly(vplot)
  })

  ## Download the displayed plot
  output$downloadPlot <- downloadHandler(
    filename = function() { paste("currentPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, width=10, height=6, plot = vplot, device = "pdf")
    }
  )

  #######################################
  ## Search for co-eluting proteins tab #
  #######################################

  output$baseProt <- renderUI({
    selectInput("baseProtein", "Base Protein", input$fvalue)
  })

  observeEvent(input$plot1_brush, {
    apex$bound_left <- round(input$plot1_brush$xmin,digits = 0)
    apex$bound_right <- round(input$plot1_brush$xmax,digits = 0)
  })

  ## Plots for semi-targeted search
  # Selection plot

  output$plot1 <- renderPlot({
    ## repl <- as.numeric(gsub(".*_r","",input$trace))
   ## cond <- ifelse(gsub("_r.*","",input$trace) == "mit", "Mitosis", "Interphase")
    p <- ggplot(target_id_traces(), aes(x=fraction, y=intensity_mean)) + #[replicate == repl]
      xlab("SEC fraction number(apparent MW[kDa])")+
      ylab("Protein level SWATH-MS intensity (top2 peptide sum)") +
      geom_line(aes(color = Gene_names, linetype = condition), size = 1, alpha = 0.8) +
      theme_bw() +
      ggtitle(unique(target_id_traces()$Gene_names)) +
      scale_x_continuous(breaks = lx.frc, labels = lx)

    if (!is.null(apex$bound_right)){
      p <- p + geom_vline(xintercept=apex$bound_left,linetype="dashed")
      p <- p + geom_vline(xintercept=apex$bound_right,linetype="dashed")
    }

    if(input$show_monomers){
      p <- p + geom_point(aes(x = monomer_fraction, color = Gene_names, y = 0), shape = 23, fill = "white", size = 3)
    }

    p <- p +
      geom_line(data = target_id_traces()[condition == input$trace & eval(as.name(input$fcolumn)) == input$baseProtein],
                aes_string(x='fraction', y='intensity_mean', color='Gene_names', linetype = 'condition'), lwd=2)

    p
  })
  output$plots <- renderUI({
    plotOutput("plot1",
               dblclick = "plot1_click",
               brush = brushOpts(id = "plot1_brush",
                                 direction = "x",
                                 resetOnNew = TRUE))
  })

  observeEvent(input$reset, {
    output$plots <- renderUI({
      plotOutput("plot1",
                 dblclick = "plot1_click",
                 brush = brushOpts(id = "plot1_brush",
                                   direction = "x",
                                   resetOnNew = TRUE))

    })
  })

  # Result plot

  observeEvent(input$search, {
    # Make download button appear
    output$downloadData <- downloadHandler(
      filename = function(){ paste0("Tageted_search_result_", input$baseProtein,".csv")},
      content = function(file){
        write.csv(file = file, x = searchResFilt(), quote = F, row.names = F)
      }
    )
    output$downloadSemiTargetedRes <- renderUI({downloadButton('downloadData', 'Download')})

    # Render a new Plot
    output$plots <- renderUI({
      plotlyOutput("plot2")
    })
    output$sliderloc <- renderUI({
      sliderInput("corr", "Correllation", 0, 1, value = 0.8, step = 0.01)
    })
    output$sliderglob <- renderUI({
      sliderInput("globalCorr", "Global Correllation", 0, 1, value = 0, step = 0.01)
    })
    # Perform the search
    target_id <- up[which(up[[input$fcolumn]] == input$baseProtein),
                                      unique(Entry)]
    ## traces <- get(paste0("prot_", input$trace))
    # We take the mean intensity so the replicate does not matter
    traces <- shinySemiTargetedSearchAdapter(tr[condition == "Mitosis" & replicate == 1])
    searchRes <<- searchSemiTargeted(traces,
                       Id = target_id,
                       lower_bound = apex$bound_left,
                       upper_bound = apex$bound_right)
    searchResFilt(searchRes[cor >= 0.8 & global_cor >= 0])

    # Plot the search result
    output$plot2 <- renderPlotly({
      p <- plotSemiTargeted(search_result = searchResFilt(),
                            traces = traces,
                            Id = target_id,
                            nr_traces = 20,
                            lower_bound = apex$bound_left,
                            upper_bound = apex$bound_left,
                            plot = F)

      p <- p +
        xlab("SEC fraction number(apparent MW[kDa])")+
        ylab("Protein intensity (top2 peptide sum)") +
        scale_x_continuous(breaks = lx.frc, labels = lx)

      if (!is.null(apex$bound_right)){
        p <- p + geom_vline(xintercept=apex$bound_left,linetype="dashed")
        p <- p + geom_vline(xintercept=apex$bound_right,linetype="dashed")
      }
      ggplotly(p)
    })

  # Watch the slider input
    observeEvent({
      input$corr
      input$globalCorr}, {
        if(!is.null(input$corr)){
          searchResFilt(searchRes[cor >= input$corr])
        }
        if(!is.null(input$globalCorr)){
          searchResFilt(searchResFilt()[global_cor >= input$globalCorr])
        }
      })
    ## Table output
    output$table <- DT::renderDataTable({
      target_id <- up[which(up[[input$fcolumn]] %in% input$fvalue),
                                      unique(Entry)]
      up[Entry %in% target_id]
    })
    output$restable <- DT::renderDataTable({
      as.data.frame(searchResFilt())
    })
  })

  ## Watch the pasteInteractors button for the semitargeted search
  observeEvent(input$pastediff, {
      if(!is.null(input$corr)){
          searchResFilt(searchRes[cor >= input$corr])
      }
      if(!is.null(input$globalCorr)){
          searchResFilt(searchResFilt()[global_cor >= input$globalCorr])
      }
      paste_ids <- up[Entry %in% searchResFilt()$id][[input$fcolumn]]
      
      print(ids)

      updateSelectizeInput(session, "fvalue", selected = ids)
  })

  # Download content
  output$downloadSemiTargetedRes <- renderUI("")

  ##################################
  ## Query String interactors tab  #
  ##################################

  output$stringProt <- renderUI({
    selectInput("stringProtein", "Search Protein", input$fvalue)
  })
  output$nrinteractors <- renderUI({
    sliderInput("nrint", "Nr interactors", 0, 20, value = 10, step = 1)
  })
  output$confidencethr <- renderUI({
    sliderInput("conf", "Confidence threshold", 0, 1000, value = 400, step = 5)
  })

  ## Create the output table
  string_table <- eventReactive(target_string(),{
      restable <- stringLinks[protein1 == target_string() | protein2 == target_string()]
      restable$Protein1 <- stringIdMap[restable$protein1]
      restable$Protein2 <- stringIdMap[restable$protein2]
      restable <- restable[order(combined_score, Protein1, Protein2, decreasing = T),.(Protein1, Protein2, combined_score)]
      ## strids <- unique(restable$Protein1)[1:input$nrint] #Make the ids available globally
      restable
  })

  ## collect chromatograms of all string interactors
  string_id_traces <- eventReactive(string_table(),{
      string_ids <- unique(string_table()$Protein1)[1:input$nrint]
      ## target_ids <- unique(up[Gene_names %in% string_ids][[input$fcolumn]])
      target_id_traces <- tr[Gene_names %in% string_ids]
      target_id_traces[, monomer_fraction:=calibration_functions$MWtoFraction(protein_mw)]
      target_id_traces
  })

    
  # Plot the String interaction of selected Protein
  target_string <- eventReactive(input$stringProtein,{
    target_id <- up[which(up[[input$fcolumn]] == input$stringProtein),
                    unique(Entry)]
    ## print(target_id)
    target_string <- up[Entry == target_id]$`Cross-reference_(STRING)`
    target_string <- strsplit(target_string, split = ";")[[1]]
    ## print(target_string)
  })

  ## Render the string image
  output$plot_string_neighbors <- renderImage({

    outfile <- tempfile(fileext = '.svg')
    print(outfile)
    obtainNeighborImage(target_string(), required_score = input$conf,
                        add_white_nodes = input$nrint, type = "svg",
                        network_flavor = "evidence", filename = outfile, verbose = T)
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/svg+xml',
         # width="100%",
         height="100%",
         alt = "Cannot display string network")
  }, deleteFile = F)

  # Plot the string interactor traces
  output$stringtraces <- renderPlotly({
      vplot <- viewerPlot(string_id_traces(), split=F, errbars=F)
      ggplotly(vplot)
  })

  # String table output
  output$stringtable <- DT::renderDataTable({
    string_table()
  })

  ## Watch the pasteInteractors button for the string interactors
  observeEvent(input$paste, {
      stringids <- unique(string_id_traces()$id)
      paste_ids <- up[Entry %in% stringids][[input$fcolumn]]
      ids <- unique(c(input$fvalue, paste_ids))
      print(ids)

      updateSelectizeInput(session, "fvalue", selected = ids)
  })

  ######################################
  ## View differential Association tab #
  ######################################

  # Plot the differentially expressed proteins between conditions

  output$plot_diffexpr <- renderPlotly({
    setkeyv(diffExprProt, input$fcolumn)
    # print(trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] %in% input$fvalue)]$protein_id)
    highlight_proteins <- diffExprProt[diffExprProt[[input$fcolumn]] %in% input$fvalue]
    p <- ggplot(as.data.frame(diffExprProt), aes(x = medianLog2FC, y = -log10(pBHadj))) +
      geom_point(aes_string(group=input$fcolumn)) +
      geom_hline(yintercept = -log10(0.05)) +
      geom_vline(xintercept = c(1, -1)) +
      geom_point(data = highlight_proteins,
                 aes_string(color=input$fcolumn), size=2) +
      theme_bw()

    ggplotly(p, source = "plot_diffexpr") %>% layout(dragmode = "select")
  })


  output$volcanotraces <- renderPlotly({
    vplot <- viewerPlot(target_id_traces(), split=F, errbars=F)
    ggplotly(vplot)
  })

  observeEvent(event_data("plotly_click", source="plot_diffexpr"),{
    ed <- event_data("plotly_click", source="plot_diffexpr")
    ## print(ed$pointNumber)
    selected <- diffExprProt[[input$fcolumn]][ed$pointNumber + 1]
    ## print(selected)
    ids <- c(setdiff(input$fvalue, selected), setdiff(selected, input$fvalue))
    ## print(ids)
    updateSelectizeInput(session, "fvalue", selected = ids)
  })

    observeEvent(event_data("plotly_selected", source="plot_diffexpr"),{
        ed <- event_data("plotly_selected", source="plot_diffexpr")
        ## print(ed)
        ## print(ed$pointNumber)
        selected <- diffExprProt[[input$fcolumn]][ed[ed$curveNumber==0,]$pointNumber + 1]
        ## For dragging we only want to unselect if all the points in the area
        ## are selected. This is more intuitive than always having the difference
        if(length(setdiff(selected,input$fvalue))==0){
            ids <- setdiff(input$fvalue, selected)
        }else{
            ids <- unique(c(input$fvalue, selected))
        }
        ## print(ids)
        updateSelectizeInput(session, "fvalue", selected = ids)
    })

})



