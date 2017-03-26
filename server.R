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

# server definition

shinyServer(function(input, output) {
  
  ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
  output$fcolumnvalues <- renderUI({
    values <- sort(unique(get(input$fcolumn)))
    # values <- values[nchar(values)>0]
    selectizeInput("fvalue", "Search or select protein of interest", values, multiple = TRUE, options = list(maxOptions = 6000))
  })
  
  ## generate selected protein SEC traces plot
  output$plot <- renderPlotly({
    
    # collect chromatograms for selection
    ######################################
    ######################################
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] %in% input$fvalue),
                                      unique(protein_id)]
    
    
    # interphase
    target_id_traces_int_r1 <- toLongFormat(subset(prot_int_r1, trace_ids = target_id)$traces)
    target_id_traces_int_r1$replicate <- 1
    target_id_traces_int_r1$condition <- "Interphase"
    
    target_id_traces_int_r2 <- toLongFormat(subset(prot_int_r2, trace_ids = target_id)$traces)
    target_id_traces_int_r2$replicate <- 2
    target_id_traces_int_r2$condition <- "Interphase"
    
    target_id_traces_int_r3 <- toLongFormat(subset(prot_int_r3, trace_ids = target_id)$traces)
    target_id_traces_int_r3$replicate <- 3
    target_id_traces_int_r3$condition <- "Interphase"
    
    # mitosis
    target_id_traces_mit_r1 <- toLongFormat(subset(prot_mit_r1, trace_ids = target_id)$traces)
    target_id_traces_mit_r1$replicate <- 1
    target_id_traces_mit_r1$condition <- "Mitosis"
    
    target_id_traces_mit_r2 <- toLongFormat(subset(prot_mit_r2, trace_ids = target_id)$traces)
    target_id_traces_mit_r2$replicate <- 2
    target_id_traces_mit_r2$condition <- "Mitosis"
    
    target_id_traces_mit_r3 <- toLongFormat(subset(prot_mit_r3, trace_ids = target_id)$traces)
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
    
    # PLOT
    ########################
    ########################
    lx.frc <- seq(5,(ncol(prot_int_r1$traces)-1),5)
    lx <- paste(lx.frc , round(calibration_functions$SECfractionToMW(lx.frc), 1) , sep = '(' )
    lx <- paste(lx, "", sep = ')' )
    
    p <- ggplot(target_id_traces[replicate == input$replicate], aes(x=fraction, y=intensity)) +
      xlab("SEC fraction number(apparent MW[kDa])")+
      ylab("Protein level SWATH-MS intensity (top2 peptide sum)") 
    
    if (input$split_plot){
      p <- p + geom_line(aes(group = Gene_names, color = Gene_names), size = 2, alpha = 0.8) +
        theme_bw() + 
        ggtitle(unique(target_id_traces$Gene_names)) +
        scale_x_continuous(breaks = lx.frc, labels = lx)  + 
        facet_wrap(~condition, ncol = 1)
    } else{
      p <- p + geom_line(aes(color = Gene_names, linetype = condition), size = 2, alpha = 0.8) +
        theme_bw() + 
        ggtitle(unique(target_id_traces$Gene_names)) +
        scale_x_continuous(breaks = lx.frc, labels = lx) 
    }
      
    if(input$show_monomers){
      p <- p + geom_point(aes(x = monomer_fraction, color = Gene_names, y = 0), shape = 23, fill = "white", size = 3) 
    }
    
    if (input$logscale){
      p <- p + scale_y_log10()
    }
    
    ggplotly(p)
    dev.off()
    p
  })
  
  # Table output
  output$table <- renderDataTable({
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] %in% input$fvalue),
                                      unique(protein_id)]
    trace_annotation_cum[protein_id %in% target_id]
  })
  
})

