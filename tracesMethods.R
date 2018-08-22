#' Subset a traces.obj by trace_ids or fraction_ids.
#' @param traces.obj An object of type \code{traces.obj}.
#' @param trace_ids A character vector specifying the trace identifiers
#'        for subsetting \code{traces.obj}.
#' @param fraction_ids A numeric vector specifying the fraction identifiers
#'        for subsetting \code{traces.obj}.
#' @return traces.obj An object of type \code{traces.obj}.
#' @export
subset.traces <- function(traces,trace_subset_ids=NULL,trace_subset_type="id",fraction_ids=NULL){
  if (!is.null(trace_subset_ids)) {
    if (trace_subset_type %in% names(traces$trace_annotation)) {
      traces$trace_annotation <- subset(traces$trace_annotation, get(trace_subset_type) %in% trace_subset_ids)
      trace_ids <- traces$trace_annotation$id
      traces$traces <- subset(traces$traces,id %in% trace_ids)
      if (nrow(traces$traces) == 0) {
        message("Caution! Subsetting returns empty traces object.")
      }
    } else {
      stop(paste0(trace_subset_type, "is not a valid trace_subset_type."))
    }
  }
  if (!is.null(fraction_ids)){
    traces$traces <- subset(traces$traces,select=c(names(traces$traces)[fraction_ids],"id"))
    traces$fraction_annotation <- subset(traces$fraction_annotation ,id %in% fraction_ids)
    if (nrow(traces$traces) == 0) {
      stop(paste0("fraction_ids (",fraction_ids,") do not match the available fractions in the traces object."))
    }
  }
  traces
}

#' Get a matrix of intensity values from a traces object.
#' @param traces.obj An object of type \code{traces.obj}.
#' @return A matrix with intensity values.
#' @export
getIntensityMatrix <- function(traces.obj) {
    ids <- traces.obj$traces$id
    intensity.mat <- as.matrix(sapply(subset(traces.obj$traces,
                                      select=-id),as.numeric))
    rownames(intensity.mat) <- ids
    intensity.mat
}


#' Convert a data.table containing traces from wide format to long format.
#' @param traces.dt A data.table with an id column \code{id} and
#'        columns of continuously numbered fractions.
#' @return A data.table with columns
#'          \itemize{
#'           \item \code{id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
toLongFormat <- function(traces.dt) {
  traces.dt.long <-
    melt(traces.dt, id.var='id', variable.name='fraction',
         value.name='intensity', variable.factor=FALSE)
  traces.dt.long[, fraction := as.numeric(fraction)]
  setkey(traces.dt.long,id)
  data.table(traces.dt.long)
  traces.dt.long
}

#' Plot a traces.obj.
#' @param traces.obj An object of type \code{traces.obj}.
#' @export
plot.traces <- function(traces.obj, plot=TRUE, ledgend = TRUE, title=NULL) {
  traces.long <- toLongFormat(traces.obj$traces)
  pl <- ggplot(traces.long)
  if (is.null(title)) {
    pl <- pl + ggtitle(paste(traces.obj$trace_type, 'traces'))
  } else {
    pl <- pl + ggtitle(title)
  }
  pl <- pl + xlab('fraction') + ylab('intensity')
  pl <- pl + geom_line(aes(x=fraction, y=intensity, color=id))
  if (!ledgend) {
    pl <- pl + theme(legend.position="none")
  }
  if (plot) plot(pl)
  pl
}


#' Summarize a \code{traces.obj}
#' @param traces.obj An object of type \code{traces.obj}.
#' @export
summary.traces <- function(traces.obj) {
  no_traces <- nrow(traces.obj$trace_annotation)
  no_decoys <- length(grep("DECOY", traces.obj$trace_annotation$protein_id))
  no_targets <- no_traces - no_decoys
  pct_decoys <- signif(no_decoys/no_traces * 100, 2)
  res <- c(no_traces, no_targets, no_decoys, pct_decoys)
  names(res) <- c("No. of Traces", "No. of Targets", "No. of Decoys", "% Decoys")
  res
}

#' Test if an object is of class traces.
#' @param traces Object of class traces.
#' @param type Character string specifying whether a specific type of traces is required.
#' The two options are "peptide" or "protein". Default is \code{NULL},
#' meaning that no specific type is required.
.tracesTest <- function(traces,type=NULL){
  if (! class(traces)=="traces") {
    stop("Object is not of class traces.")
  }
  if (! all(names(traces)==c("traces","trace_type","trace_annotation","fraction_annotation"))) {
    stop("Traces object doesn't contain all necessary items: traces, trace_type, trace_annotation, and fraction_annotation.")
  }
  if (!is.null(type)) {
    if (type != traces$trace_type) {
      stop("Traces object is of wrong type. Please check your input traces.")
    }
  }
  if (! identical(traces$traces$id,traces$trace_annotation$id)) {
    stop("IDs in traces and trace_annotation are not identical.")
  }
  if (! identical(names(traces$traces),c(traces$fraction_annotation$id,"id"))) {
    stop("Fractions in traces and fraction_annotation are not identical.")
  }
}

