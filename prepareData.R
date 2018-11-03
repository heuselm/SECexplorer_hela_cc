# prepare Data for visualization in shiny web app
# This script can be run from the folder where the shiny program is saved and
# will create the correct data folder

# Load libraries
library(data.table)
library(devtools)
devtools::install_github("CCprofiler/CCprofiler", ref = "helaCC")
library(CCprofiler)

## Setup variables
shinydir <- "~/code/SECexplorer_hela_cc/"

############################
## Generate the app data
############################

## Setup the calibration of SEC fraction to MW
calibrationTable <- data.table(std_weights_kDa = c(1398, 699, 300, 150, 44, 17),
                               std_elu_fractions = c(15.25, 24, 31, 39.275, 47.5, 53.5)
                               )

calibration_functions = calibrateMW(calibration_table = calibrationTable,
                                    plot=FALSE,
                                    PDF=FALSE)

## Import the protein traces
# Obtained with script: 01_proteinInferenceAndPlotting.R
protTraces <- readRDS("~/home/HeLa_CCprofilerAnalysis/results/proteinInferenceAndPlotting/proteinTracesLong_mean_sd_sem.rda")

## Download annotation data from uniprot
# This is an annotation table of the human swissprot assembly 
up <- fread("~/sonas/databases/Uniprot/uniprot-all-human9606-20161130_extended.tab")
up[, Mass:=as.numeric(gsub(",", ".", Mass))]
names(up) <- gsub(" ", "_", names(up))
up <- up[Entry %in% protTraces$id]

# Select the id columns that the user can use to select a protein
idcols <- names(up)[c(1, 4, 5, 13, 19:23)]

## Create a password (Do not put this section with the actual pw under public version control!)

pwd <- "mypassword"

############################
## Create the data directory
############################

setwd(shinydir)
dir.create("www/data")
setwd("www/data")

saveRDS(calibration_functions, "calibration_functions.rda")
saveRDS(up, "uniprotMapping.rda")
saveRDS(protTraces, "proteinTracesLong_mean_sd_sem.rda")
saveRDS(pwd, "pass.rda")
saveRDS(idcols, "idcols.rda")
