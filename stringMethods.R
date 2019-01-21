obtainNeighborImage <- function (identifiers, required_score = 400, add_white_nodes = 10, type = c("png","svg"),
                         network_flavor = "evidence", filename = NULL, str_version = 10, verbose = F){
  type <- match.arg(type)
  identifiers <- unique(identifiers)
  identifiers <- identifiers[!is.na(identifiers)]
  identifiers <- identifiers[identifiers != ""]
  if (length(identifiers) > 400) {
    stop("Input is >400 id's long")
  }

  urlStr <- paste0("http://string-db.org/api/", type,"/network")

  identifier_str <- paste(identifiers, collapse = "%0D")

  urlStr <- paste0(urlStr, "?identifiers=", identifier_str)

  # Paste the parameters
  urlStr <- paste0(urlStr, "&")
  urlStr <- paste0(urlStr, "add_white_nodes=", add_white_nodes)
  urlStr <- paste0(urlStr, "&")
  urlStr <- paste0(urlStr, "network_flavor=", network_flavor)
  urlStr <- paste0(urlStr, "&")
  urlStr <- paste0(urlStr, "required_score=", required_score)
  if(verbose) message(urlStr)
  if (!is.null(file)){
    download.file(urlStr, destfile = filename, quiet = !verbose)
    # img <- readPNG(filename)
  } else {
    download.file(urlStr, destfile = "tmp.png", quiet = !verbose)
    # img <- readPNG("tmp.png")
  }
  # return(img)
}
