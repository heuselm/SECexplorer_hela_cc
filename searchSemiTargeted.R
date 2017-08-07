#' 
#' Search for local co-elution of traces
#' @description Searches for highly correlated protein elution profiles,
#' based on a protein/peptide trace of choice
#' @import data.table
#' @param traces Object of class traces.
#' @param Id Character string, the id of the protein/peptide to be searched
#' @param lower_bound Numeric, the lowest fraction number of the region to be correlated.
#' If \code{NULL}, the lowest possible fraction is taken.
#' @param upper_bound Numeric, the highest fraction number of the region to be correlated
#' If \code{NULL}, the highest possible fraction is taken.
#' @param corr_method Character string indicating which correlation coefficient is to be used.
#' One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param plot Logical, wether to plot a visualization of the search result, defaults to 
#' \code{FALSE}.
#' @param nr_traces Numeric integer, the top x traces from the semi-targeted search
#' which should be plotted. (Only applied if plot is \code{TRUE}, defaults to 10)
#' @return data.table with interaction candidates ranked based on their correllation.
#' Has the following columns:
#' \itemize{
#' \item 'id' holds trace Ids
#' \item 'cor' holds the local correllation in the range specufied by the boundaries.
#' \item 'global_cor' holds the correlation over all fractions
#' \item Additional information columns from the trace annotation
#' }
#' @export
#' @examples 
#' ## Load example Data
#' proteinTraces <- exampleProteinTraces
#' 
#' ## Perform the semi-targeted search for protein P12956 between fractions 38 trough 47
#' searchList <- searchSemiTargeted(traces = proteinTraces,
#'                                  Id = "P12956",
#'                                  lower_bound = 38,
#'                                  upper_bound = 47)
#' ##Inspect the result
#' searchList
#' ## The 'cor' colum holds the local correllation in the range specufied by the boundaries.
#' ## The 'global_cor' column holds the correlation over all fractions

searchSemiTargeted <- function(traces,
                               Id,
                               lower_bound = NULL,
                               upper_bound = NULL,
                               corr_method = "pearson",
                               plot = FALSE,
                               nr_traces = 10){
  
  if(is.null(lower_bound)) lower_bound <- min(traces$fraction_annotation$id)
  if(is.null(upper_bound)) upper_bound <- max(traces$fraction_annotation$id)
  if((!is.numeric(lower_bound)) | !(lower_bound %in% traces$fraction_annotation$id)){
    stop("lower_bound must be equal to an existing fraction")
  }
  if((!is.numeric(upper_bound)) | !(upper_bound %in% traces$fraction_annotation$id)){
    stop("upper_bound must be equal to an existing fraction")
  }
  if(!any(traces$trace_annotation$id == Id)){
    if(traces$trace_type == "peptide"){
      if(!any(traces$trace_annotation$protein_id == Id)){
        stop(paste("Id: ", Id, "not found."))
      }else{
        ids <- traces$trace_annotation[protein_id == Id, id]
        res <- lapply(ids, searchSemiTargeted, traces = traces, lower_bound = lower_bound,
                      upper_bound = upper_bound, corr_method = corr_method, plot = FALSE)
        for(i in 1:length(res)){
          res[[i]]$base_protein <- Id
          res[[i]]$base_peptide <- ids[i]
        }
        resTable <- do.call('rbind', res)
        return(resTable)
      }
    }else{
      stop(paste("Id: ", Id, "not found."))
    }
  } 
  
  quantm <- getIntensityMatrix(traces = traces)
  search_trace <- quantm[Id,]
  quantm <- quantm[!(rowSums(quantm) == 0),]
  global_cor <- cor(t(quantm), search_trace, method = corr_method)
  quantm <- quantm[, lower_bound:upper_bound]
  search_trace <- search_trace[lower_bound:upper_bound]
  res <- cor(t(quantm), search_trace, method = corr_method)
  res <- data.table(id = rownames(res),
                    cor = res[,1],
                    global_cor = global_cor[,1])
  res_info <- merge(res, traces$trace_annotation, by = "id", all.x = T)
  res_info <- res_info[order(cor, decreasing = TRUE)]
  # setnames(res_info, "cor", paste0("cor", lower_bound, "_", upper_bound))
  if(plot){
    plotSemiTargeted(search_result = res_info,
                     traces = traces,
                     Id = Id,
                     nr_traces = nr_traces,
                     lower_bound = lower_bound,
                     upper_bound = upper_bound)
  }
  return(res_info)
}



#' 
#' Plot locally co-eluting traces
#' @description Plot locally co-eluting traces resulting from a semi-targeted search
#' @import data.table
#' @import ggplot2
#' @param search_result a data.table with local correlations of traces as 
#' returned by \link{searchSemiTargeted}.
#' @param traces Object of class traces.
#' @param Id Character string, the id of the protein/peptide to be searched
#' @param lower_bound Numeric, the lowest fraction number of the region to be correlated.
#' If \code{NULL}, the lowest possible fraction is taken.
#' @param upper_bound Numeric, the highest fraction number of the region to be correlated
#' If \code{NULL}, the highest possible fraction is taken.
#' @param nr_traces Numeric integer, the top x traces from the semi-targeted search
#' which should be plotted.
#' @return A ggplot of locally corellating traces.
#' @export
#' @examples 
#' ## Load example Data
#' proteinTraces <- exampleProteinTraces
#' 
#' ## Perform the semi-targeted search for protein P12956 between fractions 38 trough 47
#' searchList <- searchSemiTargeted(traces = proteinTraces,
#'                                  Id = "P12956",
#'                                  lower_bound = 38,
#'                                  upper_bound = 47)
#' ## Plot the result
#' plotSemiTargeted(search_result = searchList,
#'                  traces = proteinTraces,
#'                  Id = "P12956",
#'                  nr_traces = 10,
#'                  lower_bound = 38,
#'                  upper_bound = 47)
#' 
#' 
plotSemiTargeted <- function(search_result,
                             traces,
                             Id,
                             nr_traces = 10,
                             lower_bound = NULL,
                             upper_bound = NULL,
                             plot = T,
                             PDF = FALSE,
                             name = "semiTargetedSearchPlot"){
  
  # if (traces$trace_type == "protein") {
  #   feature <- data.table(complex_id = id,
  #                         subunits_with_signal = paste(search_result[1:nr_traces, id], sep = ";"),
  #                         )
  # }else if(traces$trace_type == "peptide"){
  #   feature <- data.table(protein_id = traces$trace_annotation[id == id, protein_id])
  #   
  # }else{
  #   stop("Traces object must be of type 'peptide' or 'protein'")
  # }
  subset_ids <- unique(search_result$protein_id)
  nr_hits <- length(subset_ids)
  subset_ids <- subset_ids[1:min(nr_traces, nr_hits)]
  tracesSubs <- subset(traces, trace_subset_ids = subset_ids, trace_subset_type = "protein_id")
  tracesLong <- toLongFormat(tracesSubs$traces)
  tracesLong <- merge(tracesLong, search_result[, .(id, cor, protein_id, Gene_names)], all.x = T)
  
  p <- ggplot(tracesLong) +
    geom_line(data = tracesLong[id != Id],
              aes_string(x='fraction', y='intensity',
                         group = 'id', color='Gene_names', size = 'cor')) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    geom_line(data = tracesLong[id == Id],
              aes_string(x='fraction', y='intensity', group = 'id'), lty = "dashed") +
    scale_size(range = c(0,2))
  if(!is.null(upper_bound) & !is.null(lower_bound)){
    p <- p +
      geom_vline(xintercept=lower_bound, colour="grey",linetype="longdash", alpha=0.4) +
      geom_vline(xintercept=upper_bound, colour="grey",linetype="longdash", alpha=0.4) +
      geom_rect(xmin = lower_bound, xmax = upper_bound, ymin = -Inf, ymax = Inf, fill="grey", alpha = 0.15)
  }
  if(plot){
    if(PDF){
      pdf(gsub("\\.pdf$|$",".pdf", name))
    }
    plot(p)
    if(PDF){
      dev.off()
    }
  }else{
    return(p)
  }
}

#' 
#' Filter a semi-targeted search result
#' @description Uses the sibling peptide correllation of a protein of interest to
#' filter proteins that locally co-elute with that protein.
#' @import data.table
#' @import ggplot2
#' @param search_result a data.table with local correlations of traces as 
#' returned by \link{searchSemiTargeted}. Must come from a semi-targeted search that
#' was performed on peptideTraces with a protein Id.
#' @param Id Character string, the id of the protein that was searched
#' @param zscore_cutoff Numeric, Specifies how many standard deviations below the mean
#' sibling peptide correlation the cutoff for filtering is to be set.
#' @details The mean correlation between all peptides of the target protein is
#' calculated (mean sibling peptide correllation) as well as the mean correllations
#' of all peptides of the target protein with peptides of all other proteins. 
#' @return A ggplot of locally corellating traces.
#' @export

filterSemiSearchResult <- function(search_result,
                                   Id,
                                   zscore_cutoff = 1,
                                   returnProtLevel = F){

  sm <- search_result[, .(mean = mean(cor, na.rm = T),
                          sd = sd(cor, na.rm = T),
                          median = median(cor, na.rm = T),
                          globalMean = mean(cor, na.rm = T),
                          globalSd = sd(cor, na.rm = T),
                          GlobalMedian = median(cor, na.rm = T)),
                      by = protein_id]
  cutoff = sm[protein_id == Id, mean - zscore_cutoff * sd]
  
  smFilt <- sm[mean >= cutoff]
  res <- search_result[protein_id %in% smFilt$protein_id]
  return(res)
}



