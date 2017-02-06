source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
library(GEOquery)
Set working directory for download
setwd("/Users/ogriffit/Dropbox/BioStars")

#Download the CEL file package for this dataset (by GSE - Geo series id)
getGEOSuppFiles("GSE27447")

#Unpack the CEL files
setwd("/Users/ogriffit/Dropbox/BioStars/GSE27447")
untar("GSE27447_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")
