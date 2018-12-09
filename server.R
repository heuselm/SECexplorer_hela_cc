library(shiny)
library(ggplot2)
library(plotly)
library(data.table)
library(STRINGdb)
canCrop <- require(magick) # This package is not available on windows (used to crop images)


# data preparation
setwd("www/data")
pass <- readRDS("pass.rda")
## trace_annotation_cum <- readRDS("trace_ann.rda")
calibration_functions <- readRDS("calibration_functions.rda")
up <- readRDS("uniprotMapping.rda")
## load("data_.rda")
trall <- readRDS("proteinTracesLong_mean_sd_sem.rda")
tr <- trall[id %in% c("P37198", "Q7Z3B4", "Q9BVL2")]
stringLinks <- fread("9606.protein.links.v10.5.HeLaSubset.txt")
stringIdMap <- readRDS("stringIdMapUniq.rda")
#Load the differential expression data
# diffExprProt <- readRDS("differentiallyExpressedProteins.rda")
# diffExprProt <- merge(diffExprProt, trace_annotation_cum, by.x = "feature_id", by.y = "protein_id")
# saveRDS(diffExprProt, "differentiallyExpressedProteinsAnn.rda")
diffExprProt <- readRDS("differentiallyExpressedProteinsAnn.rda")
setwd("../../")
default_proteins <- c("NUP54", "NUP62", "NUP58 KIAA0410 NUPL1")

# Source the functions
source("searchSemiTargeted.R")
source("tracesMethods.R")
source("stringMethods.R")
# string_db <- STRINGdb$new( version="10", species=9606, score_threshold=400, input_directory="" )

idcols <- readRDS("www/data/idcols.rda")
## for (i in seq_along(idcols)){
##   assign(idcols[i], up[idcols][[i]])
## }

# server definition

shinyServer(function(input, output, session) {

  ## Set some global variables
  apex <- reactiveValues(bound_left = NULL, bound_right = NULL)
  searchResFilt <- reactiveVal(NULL)

  ############################
  ## Viewer Tab
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
                     choices=c("Standard Error (SEM)", "Standard deviation (SD)"),
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
  output$plot <- renderPlotly({
    # PLOT
    p <- ggplot(target_id_traces(), aes(x=fraction, y=intensity_mean + 1)) +
      xlab("SEC fraction number(apparent MW[kDa])") +
      ylab("Protein level SWATH-MS intensity (top2 peptide sum)")

    if (input$split_plot){
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
    if(input$error_bars){
      if(input$ertype == "Standard deviation (SD)"){
        p <- p + geom_point(aes(color = Gene_names)) +
          geom_errorbar(aes(ymin=intensity_mean-intensity_sd,
                            ymax = intensity_mean+intensity_sd,
                            color = Gene_names),
                        position = position_dodge(0.2))
      }
      if(input$ertype == "Standard Error (SEM)"){
        p <- p + geom_point(aes(color = Gene_names)) +
          geom_errorbar(aes(ymin=intensity_mean-intensity_se,
                            ymax = intensity_mean+intensity_se,
                            color = Gene_names),
                        position = position_dodge(0.2))
      }
    }
    if(input$show_monomers){
      p <- p + geom_point(aes(x = monomer_fraction, color = Gene_names, y = 0), shape = 23, fill = "white", size = 3)
    }

    if (input$logscale){
      p <- p + scale_y_log10()
    }

    ggplotly(p)
  })
  ###########################
  ## Search for co-eluting proteins tab
  ###########################

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

  # Watch the pasteInteractors button
  observeEvent(input$paste, {
      if(!is.null(input$corr)){
        searchResFilt(searchRes[cor >= input$corr])
      }
      if(!is.null(input$globalCorr)){
        searchResFilt(searchResFilt()[global_cor >= input$globalCorr])
     }
    # ids <- unique(c(restable$Protein1, restable$Protein2))
    ids <- unique(c(input$fcolumn, strids))
    print(ids)

    updateSelectizeInput(session, "fvalue", selected = ids)
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

  # Table output
  output$table <- renderDataTable({
    target_id <- up[which(up[[input$fcolumn]] %in% input$fvalue),
                                      unique(Entry)]
    up[Entry %in% target_id]
  })
  output$restable <- renderDataTable({
    searchResFilt()
  })
    # Download content
    output$downloadData <- downloadHandler(
      filename = function(){ paste0("Tageted_search_result_", input$baseProtein,".csv")},
      content = function(file){
        write.csv(file = file, x = searchResFilt(), quote = F, row.names = F)
      }
    )
  })
  ###########################
  ## Query String interactors tab
  ###########################

  output$stringProt <- renderUI({
    selectInput("stringProtein", "Search Protein", input$fvalue)
  })
  output$nrinteractors <- renderUI({
    sliderInput("nrint", "Nr interactors", 0, 20, value = 10, step = 1)
  })
  output$confidencethr <- renderUI({
    sliderInput("conf", "Confidence threshold", 0, 1000, value = 400, step = 5)
  })

    # Plot the String interaction networks
    # output$plot_st_string <- renderPlot({
    #   # string_ids <- string_db$map(searchResFilt(), my_data_frame_id_col_names = "id")
    #   string_ids <- up[Entry %in% searchResFilt()$id]$`Cross-reference_(STRING)`
    #   string_ids <- string_ids[string_ids != ""]
    #   print(string_ids)
    #   p <- string_db$plot_network(string_ids, add_link = T, add_summary = T)
    #   p
    # })

  # Plot the String interaction of selected Protein
  output$plot_string_neighbors <- renderImage({
    target_id <- up[which(up[[input$fcolumn]] == input$stringProtein),
                                      unique(Entry)]
    print(target_id)
    target_string <- up[Entry == target_id]$`Cross-reference_(STRING)`
    target_string <- strsplit(target_string, split = ";")[[1]]
    print(target_string)
    # target_neighbors <- string_db$get_neighbors(target_string)
    # print(target_neighbors)
    # p <- string_db$plot_network(c(target_string, target_neighbors), add_link = F, add_summary = F)
    # p
    # width  <- session$clientData$output_plot_width
    # height <- session$clientData$output_plot_height
    # mysvgwidth <- width/96
    # mysvgheight <- height/96

    outfile <- tempfile(fileext = '.svg')
    print(outfile)
    obtainNeighborImage(target_string, required_score = input$conf,
                        add_white_nodes = input$nrint, type = "svg",
                        network_flavor = "evidence", filename = outfile, verbose = T)
    # if(canCrop){
    #   tmpImg <-image_read_svg(outfile)
    #   width <- image_info(tmpImg)$width
    #   height <- image_info(tmpImg)$height
    #   tmpImgC <- image_crop(tmpImg, geometry_area(height*1.25, height, (width-height*1.25)/2,0))
    #   outfileC <- tempfile(fileext = '.png')
    #   image_write(tmpImgC, path=outfileC, format = "png" )
    #   print(outfileC)
    #   # Return a list containing the filename
    #   list(src = outfileC,
    #        contentType = 'image/png',
    #        # width="100%",
    #        height="100%",
    #        alt = "Cannot display string network")
    # }else{
    # }
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/svg+xml',
         # width="100%",
         height="100%",
         alt = "Cannot display string network")
  }, deleteFile = F)


  #####################################
  ## View differential Association tab
  ###################################

  # Plot the differentially expressed proteins between conditions

  output$plot_diffexpr <- renderPlotly({
    setkeyv(diffExprProt, input$fcolumn)
    # print(trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] %in% input$fvalue)]$protein_id)
    highlight_proteins <- diffExprProt[feature_id %in% up[which(up[[input$fcolumn]] %in% input$fvalue)]$Entry]
    p <- ggplot(diffExprProt, aes(x = medianLog2FC, y = -log10(pBHadj))) +
      geom_point(aes(group=feature_id)) +
      # geom_point(aes(color = (feature_id %in% "P49792"))) +
      geom_hline(yintercept = -log10(0.05)) +
      geom_vline(xintercept = c(1, -1)) +
      geom_point(data = highlight_proteins, aes(group=feature_id), color="red") +
      theme_bw()

    ggplotly(p, source = "plot_diffexpr") %>% layout(dragmode = "select")
  })

  # output$diffexprsel <- renderPrint({
  #   d <- event_data("plotly_selected", source = "plot_diffexpr")
  #   print(d$key)
  #   # d <- paste0(unique(trace_annotation_cum[protein_id %in% d][[input$fcolumn]]))
  #   if (is.null(d)) "Select events appear here (unhover to clear)" else d
  # })

  ## Watch the pasteSelection button
  observeEvent(input$pastediff, {
    sel <- event_data("plotly_selected", source = "plot_diffexpr")
    print(sel)
    ids <- unique(diffExprProt[sel[sel$curveNumber == 0,]$pointNumber + 1][[input$fcolumn]])
    print(ids)
    updateSelectizeInput(session, "fvalue", selected = ids)
  })

  # String table output
  output$stringtable <- renderDataTable({
    target_id <- up[which(up[[input$fcolumn]] == input$stringProtein),
                    unique(Entry)]
    target_string <- up[Entry == target_id]$`Cross-reference_(STRING)`
    target_string <- strsplit(target_string, split = ";")[[1]]

    restable <- stringLinks[protein1 == target_string | protein2 == target_string]
    restable$Protein1 <- stringIdMap[restable$protein1]
    restable$Protein2 <- stringIdMap[restable$protein2]
    restable <- restable[order(combined_score, Protein1, Protein2, decreasing = T),.(Protein1, Protein2, combined_score)]
    strids <<- unique(restable$Protein1)[1:input$nrint] #Make the ids available globally
    restable
  })

  #########################
  ## Password input
  #########################

  output$pwdfeedback <- renderText("Limited access. Enter password")
  observeEvent(input$enterpwd, {
      if(input$pwd == pass){
        tr <<- trall
        output$pwdfeedback <- renderText("Password correct")
      } else{
        tr <<- trall[id %in% c("P37198", "Q7Z3B4", "Q9BVL2")]
        output$pwdfeedback <- renderText("Password incorrect")
      }
  })
})



