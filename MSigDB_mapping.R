####Dwonload teh database from http://software.broadinstitute.org/gsea/msigdb/collections.jsp

library(GeneAnswers)
library(biomaRt)
##Read gmt file 
x<-scan("c3.all.v5.1.entrez.gmt",what="",sep="\n")
#x<-scan("c2.all.v5.1.entrez.gmt",what="",sep="\n")
y <- strsplit(x, "\t")
names(y) <- sapply(y, function(x) x[1])
y <- sapply(y, function(x) x[3:length(x)])
y <- lapply(y, `[`, -1)
for(i in 1:length(y)){
	y[[i]] <- y[[i]][-1]
}
save(y,file="human_motif_msig.dat")

#####All these data base are created from Human gene Identifer, to create for mouse one has to orthologmapping to following way
y_mouse <- list()
for(i in 1:length(y)){
	print(i)
	homoLL <-getHomoGeneIDs(y[[i]], species="human", speciesL="mouse",mappingMethod='biomaRt')
	y_mouse[[i]] <- as.character(homoLL)
}
names(y_mouse) <- names(y)
save(y_mouse,file="mouse_motif_msig.dat")


pathway anaysis download form concesus data base
y<-scan("CPDB_pathways_genes.tab",what="",sep="\n")
ndemo <- vector(mode="list")
for (k in 1:length(y)){
    x=strsplit(y[k],"\t")[[1]]
    ndemo[[x[1]]]<-strsplit(x[length(x)],",")[[1]]
    print(x[1])
}
idx <- grep("disease",names(ndemo))
cons2 <- ndemo[-idx]
save(ndemo,file="pathway.dat")




###############tf analysis
#############transcription factors
##########ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57–74 (2012)
#######http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/
library(RCurl)
download.file(
url = "https://raw.githubusercontent.com/slowkow/tftargets/master/data/tftargets.rda",
destfile = "tftargets.rda",
method = "curl"
)

# Load the file:
load("tftargets.rda")
#rownames(cnts.norm) <- as.character(mget(rownames(cnts.norm), org.Hs.egSYMBOL,ifnotfound=NA))
#mymat<-as.matrix(cnts.norm[!grepl('NA', rownames(cnts.norm)), ])

########download tf from mouse
x<-scan("c3.tft.v5.1.entrez.gmt",what="",sep="\n")
y <- strsplit(x, "\t")
names(y) <- sapply(y, function(x) x[1])
y <- sapply(y, function(x) x[3:length(x)])
y <- lapply(y, `[`, -1)
for(i in 1:length(y)){
    y[[i]] <- y[[i]][-1]
}
idx <- grep("UNKNOWN",names(y))
y <- y[-idx]
save(y,file="tft_human.dat")



https://www.biostars.org/p/18112/
http://iregulon.aertslab.org/collections.html#motifcolldownload
# Entrez Gene IDs.
TRED[["STAT3"]]
# Gene symbols used on the ITFP website.
ITFP[["STAT3"]]
# Entrez Gene IDs.
head(ENCODE[["STAT3"]], 100)
# Entrez Gene IDs.
Neph2012[["AG10803-DS12374"]][["STAT3"]]

https://www.bioconductor.org/help/workflows/generegulation/


