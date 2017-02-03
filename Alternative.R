
rm(list=ls())
library("spliceR")
library("humarray")
library(cummeRbund)

###bundd genome
require("BSgenome.Hsapiens.UCSC.hg38",character.only = TRUE)
#####Get annotated ORfs
#ucscCDS <- getCDS(selectedGenome="hg19", repoName="ensemble")

setwd("/home/asingh/script/script_human/workflow2_cuffdif_cuffcmp")
sqliteQuickSQL<-dbGetQuery
dbBeginTransaction<-dbBegin
options(ucscChromosomeNames=FALSE)
####path to merged gtf file
#gtfFile="/data/genomes/homo_sapiens/GRCh38_79/cuffcmp_GTF.GRCh38.79.gtf"

cuffDB<- readCufflinks(dir='.', rebuild=TRUE, gtfFile="/data/genomes/homo_sapiens/GRCh38_79/cuffcmp_GTF.GRCh38.79.gtf",genome="hg38",package="cummeRbund")
mySpliceRList <- prepareCuff(cuffDB, removeNonChanonicalChr=FALSE)

#mySpliceRList <- prepareCuff(cuffDB)
myTranscripts <- transcripts(mySpliceRList)
myExons <- exons(mySpliceRList)
conditions(mySpliceRList)
#####filter cufflink data 
#mySpliceRList <- preSpliceRFilter(mySpliceRList,filters=c("expressedIso", "isoOK", "expressedGenes", "geneOK"),expressionCutoff=1)
####
#ucscCDS <- getCDS(selectedGenome="hg19", repoName="ensemble")
#mySpliceRList <- annotatePTC(mySpliceRList, ucscCDS, Hsapiens,PTCDistance=50)
#save(mySpliceRList,file="mySpliceRList.dat")
########Analyze Orfs
#mySpliceRList <- annotatePTC(mySpliceRList, ucscCDS, Hsapiens,PTCDistance=50)
####### Analyze alternative splicing in transcripts
#mySpliceRList <- spliceR(mySpliceRList, compareTo='preTranscript', filters= 'isoOK')
mySpliceRList <- spliceR(mySpliceRList, compareTo='d00', filters=c('expressedGenes','geneOK', 'isoOK', 'expressedIso', 'isoClass'),expressionCutoff=1)
setwd("/home/asingh/script/script_human/Alt_splice_workflow_2")
save(mySpliceRList,file="output_spliceR.dat")

#####creat GTF file 
#generateGTF(mySpliceRList, filters="isoOK", filePrefix='./alternative_analysis') 
#generateGTF(mySpliceRList, filters=c("geneOK", "isoOK", "expressedGenes", "expressedIso"), scoreMethod="local")
#generateGTF(mySpliceRList,filters=c('expressedGenes','geneOK', 'isoOK', 'expressedIso', 'isoClass'),expressionCutoff=1)

################
####load out of splice R 
a<-load("output_spliceR.dat")

cuff <- readCufflinks()
replicates(cuff)

####read.cuff data
splicing_diff_data<- distValues(splicing(cuff))
sig_splicing_data<- subset(splicing_diff_data, (significant=='yes'))
dif<-subset(sig_splicing_data, q_value <0.01)
a<-mySpliceRList[[1]]
df<-ranged.to.data.frame(a)
df<-ranged.to.data.frame(a,,FALSE)
df<-ranged.to.data.frame(a,TRUE)
idx<-match(dif$gene_id,df$spliceR.gene_id)
resSig<-df[idx,]

#######Exon skipping/Inclusion (ESI)##72
spliceR.ESI<-resSig[!(is.na(resSig$spliceR.ESI) | resSig$spliceR.ESI=="0"), ]
ESI<-spliceR.ESI[c(1:45,53:55)]
write.xlsx(ESI,file="ESI_spliceR.xlsx",asTable = TRUE)
ESI<-sum(spliceR.ESI$spliceR.ESI)
#####Multi Exon skipping/Inclusion (MESI) ###12
spliceR.MESI<-resSig[!(is.na(resSig$spliceR.MESI) | resSig$spliceR.MESI=="0"), ]

MESI<-spliceR.MESI[c(1:44,47,53,58,59)]
write.xlsx(MESI,file="MESI_spliceR.xlsx",asTable = TRUE)
MESI<-sum(spliceR.MESI$spliceR.MESI)
#####Intron retention (ISI)###21
spliceR.ISI<-resSig[!(is.na(resSig$spliceR.ISI) | resSig$spliceR.ISI=="0"), ]
ISI<-spliceR.ISI[c(1:44,48,53,60,61)]
write.xlsx(ISI,file="ISI_spliceR.xlsx",asTable = TRUE)
ISI<-sum(spliceR.ISI$spliceR.ISI)
###########alternative splice site at A5
spliceR.A5<-resSig[!(is.na(resSig$spliceR.A5) | resSig$spliceR.A5=="0"), ]
A5<-spliceR.A5[c(1:44,49,53,62,63)]
write.xlsx(A5,file="A5_spliceR.xlsx",asTable = TRUE)
A5<-sum(spliceR.A5$spliceR.A5)


###########alternative splice site at A3
spliceR.A3<-resSig[!(is.na(resSig$spliceR.A3) | resSig$spliceR.A3=="0"), ]
A3<-spliceR.A3[c(1:44,50,53,64,65)]
write.xlsx(A3,file="A3_spliceR.xlsx",asTable = TRUE)
A3<-sum(spliceR.A3$spliceR.A3)

###########alternative splice site at ATSS ###154

spliceR.ATSS<-resSig[!(is.na(resSig$spliceR.ATSS) | resSig$spliceR.ATSS=="0"), ]
ATSS<-spliceR.ATSS[c(1:44,51,53,66,67)]
write.xlsx(ATSS,file="ATSS_spliceR.xlsx",asTable = TRUE)
ATSS<-sum(spliceR.ATSS$spliceR.ATSS)

###########alternative splice site at ATTS

spliceR.ATTS<-resSig[!(is.na(resSig$spliceR.ATTS) | resSig$spliceR.ATTS=="0"), ]
ATTS<-spliceR.ATTS[c(1:44,52,53,68,69)]
write.xlsx(ATTS,file="ATTS_spliceR.xlsx",asTable = TRUE)
ATTS<-sum(spliceR.ATTS$spliceR.ATTS)
#####MEE###no hits
spliceR.MEE<-resSig[!(is.na(resSig$spliceR.MEE) | resSig$spliceR.MEE=="0"), ]
write.xlsx(ATSS,file="MEE_spliceR.xlsx",asTable = TRUE)
MEE<-sum(spliceR.MEE$spliceR.MEE)


