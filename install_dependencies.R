# I used a conda environment to install the dependencies
# These are the the instructions for a manual install of the dependencies
install.packages(c("ggplot2", "plotly", "data.table"))
install.packages("magick") # Requires ImageMagick to be installed

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("STRINGdb", version = "3.8")
