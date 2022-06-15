## R Script for loading data

dropboxFolder <- "/root/Dropbox/"
folder <- "Progetto Stratega/140221_STRATEGA_v8_allsample/"
# folder <- "Progetto Stratega/140530_STRATEGA_v10/"
.currentFolder <- getwd()
setwd(paste0(dropboxFolder, folder))
data <- read.csv(dataset, header = TRUE, sep = ";")
setwd(.currentFolder)