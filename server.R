library(shiny)
library(ggplot2)
library(plotly)
library(data.table)
# data preparation
load("data_.rda")
source("searchSemiTargeted.R")
source("tracesMethods.R")

annotations <- names(trace_annotation_cum)
for (i in seq_along(annotations)){
  assign(annotations[i], trace_annotation_cum[[i]])
}

# server definition

shinyServer(function(input, output) {
  
  ## Set the variables
  apex <- reactiveValues(bound_left = NULL, bound_right = NULL)
  searchResFilt <- reactiveVal(NULL)
  ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
  output$fcolumnvalues <- renderUI({
    values <- sort(unique(get(input$fcolumn)))
    # values <- values[nchar(values)>0]
    selectizeInput("fvalue", "Search and select proteins of interest", values,
                   multiple = TRUE, options = list(maxOptions = 6000), selected = c("NUP54", "NUP62", "NUP58 KIAA0410 NUPL1"))
  })
  
  output$baseProt <- renderUI({
    selectInput("baseProtein", "Base Protein", input$fvalue)
  })
  
  observeEvent(input$plot1_brush, {
    apex$bound_left <- round(input$plot1_brush$xmin,digits = 0)
    apex$bound_right <- round(input$plot1_brush$xmax,digits = 0)
  })
  
  
  ## generate selected protein SEC traces plot
  ######################################
  ######################################
  lx.frc <- seq(5,(ncol(prot_int_r1$traces)-1),5)
  lx <- paste(lx.frc , round(calibration_functions$SECfractionToMW(lx.frc), 1) , sep = '(' )
  lx <- paste(lx, "", sep = ')' )
  
  # collect chromatograms for selection
  target_id_traces <- eventReactive(input$fvalue,{
    # apex$bound_right <- NULL
    # apex$bound_left <- NULL
    
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] %in% input$fvalue),
                         unique(protein_id)]
    # interphase
    target_id_traces_int_r1 <- toLongFormat(subset(prot_int_r1, trace_subset_ids = target_id)$traces)
    target_id_traces_int_r1$replicate <- 1
    target_id_traces_int_r1$condition <- "Interphase"
    
    target_id_traces_int_r2 <- toLongFormat(subset(prot_int_r2, trace_subset_ids = target_id)$traces)
    target_id_traces_int_r2$replicate <- 2
    target_id_traces_int_r2$condition <- "Interphase"
    
    target_id_traces_int_r3 <- toLongFormat(subset(prot_int_r3, trace_subset_ids = target_id)$traces)
    target_id_traces_int_r3$replicate <- 3
    target_id_traces_int_r3$condition <- "Interphase"
    
    # mitosis
    target_id_traces_mit_r1 <- toLongFormat(subset(prot_mit_r1, trace_subset_ids = target_id)$traces)
    target_id_traces_mit_r1$replicate <- 1
    target_id_traces_mit_r1$condition <- "Mitosis"
    
    target_id_traces_mit_r2 <- toLongFormat(subset(prot_mit_r2, trace_subset_ids = target_id)$traces)
    target_id_traces_mit_r2$replicate <- 2
    target_id_traces_mit_r2$condition <- "Mitosis"
    
    target_id_traces_mit_r3 <- toLongFormat(subset(prot_mit_r3, trace_subset_ids = target_id)$traces)
    target_id_traces_mit_r3$replicate <- 3
    target_id_traces_mit_r3$condition <- "Mitosis"
    
    target_id_traces <- rbind(target_id_traces_int_r1,
                              target_id_traces_int_r2,
                              target_id_traces_int_r3,
                              target_id_traces_mit_r1,
                              target_id_traces_mit_r2,
                              target_id_traces_mit_r3)
    
    target_id_traces <- merge(target_id_traces, trace_annotation_cum,
                              by.x = "id", by.y = "protein_id", all.y = FALSE, all.x = TRUE)
    
    target_id_traces[, monomer_fraction:=calibration_functions$MWtoSECfraction(protein_mw)]
    target_id_traces
  })
  
  
  output$plot <- renderPlotly({
    
    # PLOT
    p <- ggplot(target_id_traces()[replicate == input$replicate], aes(x=fraction, y=intensity)) +
      xlab("SEC fraction number(apparent MW[kDa])") +
      ylab("Protein level SWATH-MS intensity (top2 peptide sum)") 
    
    if (input$split_plot){
      p <- p + geom_line(aes(group = Gene_names, color = Gene_names), size = 1, alpha = 0.8) +
        theme_bw() + 
        ggtitle(unique(target_id_traces()$Gene_names)) +
        scale_x_continuous(breaks = lx.frc, labels = lx)  + 
        facet_wrap(~condition, ncol = 1)
    } else{
      p <- p + geom_line(aes(color = Gene_names, linetype = condition), size = 1, alpha = 0.8) +
        theme_bw() + 
        ggtitle(unique(target_id_traces()$Gene_names)) +
        scale_x_continuous(breaks = lx.frc, labels = lx) 
    }
      
    if(input$show_monomers){
      p <- p + geom_point(aes(x = monomer_fraction, color = Gene_names, y = 0), shape = 23, fill = "white", size = 3) 
    }
    
    if (input$logscale){
      p <- p + scale_y_log10()
    }
    
    ggplotly(p)
    # dev.off()
    # p
  })
  
  ## Plots for semi-targeted search
  # Selection plot
  
  output$plot1 <- renderPlot({
    repl <- as.numeric(gsub(".*_r","",input$trace))
    cond <- ifelse(gsub("_r.*","",input$trace) == "mit", "Mitosis", "Interphase")
    p <- ggplot(target_id_traces()[replicate == repl], aes(x=fraction, y=intensity)) +
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
      geom_line(data = target_id_traces()[replicate == repl & condition == cond & eval(as.name(input$fcolumn)) == input$baseProtein],
                aes_string(x='fraction', y='intensity', color='Gene_names', linetype = 'condition'), lwd=2)
    
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
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$baseProtein),
                                      unique(protein_id)]
    traces <- get(paste0("prot_", input$trace))
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
  })
  
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
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] %in% input$fvalue),
                                      unique(protein_id)]
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

