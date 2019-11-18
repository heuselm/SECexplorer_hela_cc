# Prepare string data for use in shiny app

## Download
# Go to https://string-db.org/cgi/download.pl?sessionId=IAq1BSw9rax8&species_text=Homo+sapiens
# Download 9606.protein.links.v10.5.txt.gz (65.9 Mb)
# Download 9606.protein.aliases.v10.5.txt.gz (12.3 Mb)
# Unzip both files

## Data preparation

stringLinks <- fread("9606.protein.links.v10.5.txt")
aliases <- fread("9606.protein.aliases.v10.5.txt")

helaProt <- readRDS("trace_annotation_cum.rda")
stringLinksSub <- fread("9606.protein.links.v10.5.HeLaSubset.txt")
  