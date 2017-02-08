library(leeBamViews)
library(EDASeq)
library("DESeq2")
library(Cairo)
  CairoFonts(
             regular="Arial:style=Regular",
             bold="Arial:style=Bold",
             italic="Arial:style=Italic",
             bolditalic="Arial:style=Bold Italic,BoldItalic",
             symbol="Symbol"
     )
setwd("/Users/asingh2/Desktop/endo_project/EDAanaysis/bam/")
#######readfastqc file 
files <-list.files(path = ".", pattern = "fastq", full.names = TRUE)
files<-files[c(1:6,9:25)]
names(files) <- gsub("\\.fastq.*", "", basename(files))
met <- DataFrame(conditions=c(rep("E2",2),rep("E4",2),rep("E6",2),rep("L2",2),rep("L3",2),rep("L5",2),rep("P2",2),rep("P3",2),rep("P5",2)),row.names=names(files))
fastq <- FastqFileList(files)
elementMetadata(fastq) <- met
fastq

########load bam file
files <- list.files(path=".",pattern = "bam", full.names = TRUE)
names(files) <- gsub("\\.bam", "", basename(files))
gt <- gsub(".*/", "", files)
gt <- gsub("_.*", "", gt)
lane <- gsub(".*(.)$", "\\1", gt)
geno <- gsub(".$", "", gt)
pd <- DataFrame(geno=geno, lane=lane,paste(geno,lane,sep="."))
bfs <- BamFileList(files)
elementMetadata(bfs) <- pd
bfs
colors <- c(rep(rgb(1,0,0,alpha=0.7),3),
rep(rgb(0,0,1,alpha=0.7),3),
rep(rgb(0,1,0,alpha=0.7),3),
rep(rgb(0,1,1,alpha=0.7),3))
CairoPDF(file="Aligned_reads.pdf",width=15/1.54,height=15/1.54)
barplot(bfs,las=2,col=colors,horiz=TRUE,main="Number of mapped reads")
legend("topright",unique(elementMetadata(bfs)[,1]), fill=unique(colors))
dev.off()

CairoPDF(file="Meanper-basequalityofmappedreads.pdf",width=15/1.54,height=15/1.54)
plotQuality(bfs,col=colors,lty=1)
legend("topright",unique(elementMetadata(bfs)[,1]), fill=unique(colors))
dev.off()
barplot(bfs[[1]],las=2)
dev.copy2pdf(file="Number of mapped reads per-chromosome.pdf")
plotNtFrequency(bfs[[1]])
dev.copy2pdf(file="frequencey.pdf")
###loading ocuntfile, its better to do individual 
setwd("/Users/asingh2/Desktop/endo_project/counts/")
meta<-read.delim("meta_file_all.txt",header=T,sep="")
meta<-read.delim("meta_fileLvsP.txt",header=T,sep="")
meta<-read.delim("meta_l_e.txt",header=T,sep="")
meta<-read.delim("meta_Evs_P.txt",header=T,sep="")
group_names <- meta$condition
options(stringsAsFactors = FALSE)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = meta, directory = ".",design = ~ condition)
###chekcing what is there in the DESeq2Table then type follwoing command
rowRanges(DESeq2Table)
colData(DESeq2Table)
GeneCounts <- counts(DESeq2Table)
###### Download GC content from "Ensembl biomart", attribute Ensembl gene id %GC contenet and gene length
a<-read.delim("mart_export.txt",sep=",",header=T)
tmp <- match(rownames(GeneCounts), a$Ensembl.Gene.ID)
b<-a[tmp,]
b$gc<-b$X..GC.content/100
c<-b[c(3,2,4)]
c<-na.omit(c)
###ignore the match number from the matrix 
rownames(c) <- NULL
colnames(c) <- c("Ensembl_Gene_ID","length","gc")
write.table(c, file="Mouse_GC_lengths.tsv", sep="\t")
write.csv(c, file="Mouse_GC_lengths.csv")
####Just ignor the enseml_gene_id for the next analysis 
m<-as.matrix(c[,-1])
row.names(m)<-c[,1]
m<-data.frame(m)
filter <- apply(GeneCounts,1,function(x) mean(x)>10)
table(filter)
common <- intersect(rownames(m),rownames(GeneCounts[filter,]))
setwd("/Users/asingh2/Desktop/endo_project/EDAanaysis/bam/")
common <- intersect(rownames(m),rownames(GeneCounts[filter,]))
length(common)

data <- newSeqExpressionSet(counts=as.matrix(GeneCounts[common,]),featureData=m[common,],phenoData=AnnotatedDataFrame(data.frame(conditions=factor(c(rep("L",3),rep("H",3),rep("E",3))),row.names=colnames(GeneCounts[common,]))))
head(counts(data))
pData(data)
data
#data <- newSeqExpressionSet(counts=as.matrix(GeneCounts[common,]),featureData=m[common,],phenoData=AnnotatedDataFrame(data.frame(conditions=factor(c(rep("L",3),rep("E",3))),row.names=colnames(GeneCounts[common,]))))
#head(exprs(data))
#pData(data)
#data
#data <- newSeqExpressionSet(counts=as.matrix(GeneCounts[common,]),featureData=m[common,],phenoData=AnnotatedDataFrame(data.frame(conditions=c(rep("E",3),rep("H",3)),row.names=colnames(GeneCounts[common,]))))
#head(exprs(data))
#pData(data)
#data
head(fData(data))
library(RColorBrewer)
#colors.lane <- brewer.pal(length(pData(data), "Accent")[pData(data)]
#colors.plex <- c("darkgray", "red", "blue", "orange", "lightgray")[plex]
#boxplot(data, col=colors.lane, las=2)
#boxplot(data, col=colors.plex, las=2)



colors <- c(rep(rgb(1,0,0,alpha=0.7),3),
rep(rgb(0,0,1,alpha=0.7),3),
rep(rgb(0,1,0,alpha=0.7),3),
rep(rgb(0,1,1,alpha=0.7),3))
CairoPDF(file="beforenorm.pdf",width=12/1.54,height=10/1.54)
boxplot(data,col=colors)
dev.off()
CairoPDF(file="meanvsvariance.pdf",width=12/1.54,height=10/1.54)
meanVarPlot(data, log=TRUE)
dev.off()
dataWithin <- withinLaneNormalization(data,"gc", which="full")
dataNorm <- betweenLaneNormalization(dataWithin, which="full")

CairoPDF(file="afternorm.pdf",width=12/1.54,height=10/1.54)
boxplot(dataNorm,col=colors)
dev.off()
CairoPDF(file="beforebiasgcplot.pdf",width=12/1.54,height=10/1.54)
biasPlot(data, "gc", log=TRUE, ylim=c(0, 8))
dev.off()

CairoPDF(file="afterbiasgcplot.pdf",width=12/1.54,height=10/1.54)
biasPlot(dataNorm, "gc", log=TRUE, ylim=c(0, 8))
dev.off()


#dataOffset <- withinLaneNormalization(data,"gc",which="full",offset=TRUE)
#dataOffset <- betweenLaneNormalization(dataOffset,which="full",offset=TRUE)


library(DESeq)
library(edgeR)
counts <- as(dataNorm, "CountDataSet")
counts(dataOffset)
design <- model.matrix(~ conditions, data=pData(dataOffset))
disp <- estimateGLMCommonDisp(counts(dataOffset),design, offset=-offst(dataOffset))
fit <- glmFit(counts(dataOffset), design, disp, offset=-offst(dataOffset))
lrt <- glmLRT(fit, coef=2)
summary(de <- decideTestsDGE(lrt, p=0.05, adjust="BH"))
edgeglm <- as.data.frame(topTags(lrt, n=length(rownames(fit))))
d<-subset(edgeglm, FDR <0.01)
k<-subset(d,logFC >0)
k<-subset(d,logFC <0)



out <- topTags(lrt, n=Inf, adjust.method="BH")
summary(de <- decideTestsDGE(lrt, p=0.01, adjust="BH"))

#####extract all de genes
 results <- topTags( de, n = nrow( de$table ) , sort.by = "logFC" )$table


##extract DE from edgeR 
isDE <- as.logical(de)
DEnames <- rownames(counts(dataOffset))[isDE]
geneList<-which(de[,1] == 1)
up<-lrt[geneList,]
keep <- out$table$FDR <= 0.05
keep <- out$table$FDR <= 0.05 & abs(out$table$logFC) >= 1
out[keep,]
write.table(out$table[keep,], file="DE.tsv", row.names = TRUE, col.names = NA, quote=FALSE, sep = "\t")


















#tmp <- match(rownames(GeneCounts), a$Ensembl.Gene.ID)
#merge(GeneCounts,a , by.x=c(""), by.y=c("AuthorFirstName"),all.x=TRUE)
#all <- intersect(rownames(l), rownames(k))














############
library(BSgenome)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
library(Biostrings)
library(ShortRead)
library(BSgenome.Mmusculus.UCSC.mm10)
gtf_file_path <-"/Users/asingh2/Desktop/endo_project/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"
FASTAfile="/Users/asingh2/Desktop/endo_project/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa" 
gtf <- import.gff (gtf_file_path, format = "gtf",  asRangedData = F, feature.type="exon")
grl <- reduce(split(gtf, elementMetadata(gtf)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <-  alphabetFrequency(getSeq(Mmusculus, reducedGTF), letters="CG")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")
write.table(output, file="GC_lengths.tsv", sep="\t")




data <- newSeqExpressionSet(counts=as.matrix(GeneCounts[common,]),
featureData=feature[common,],
phenoData=data.frame(
conditions=c(rep("E",3),rep("L",3),,rep("H",3)),
row.names=colnames(GeneCounts)))
##############################
#setwd("/Users/asingh2/Desktop/endo_project/EDAanaysis/bam/")
#t<-read.delim("GC_lengths.tsv",header=T,sep="")
#k<-read.delim("Mus_musculus.GRCm38.68_gene_GC.txt",header=T,sep="")
#l<-t[rownames(t) %in%k$geneid,]
#l1<-k[rownames(t) %in%k$geneid,]
#lo<-match(l1$geneid, rownames(l))
#ll<-t[lo,]
#lll<-cbind(l1,ll)
#write.table(lll, file="final_GC_lengths.tsv", sep="\t")
#a<-read.delim("final_GC_lengths.tsv",header=T,sep="")
####
#a<-read.delim("last_final_GC_lengths.csv",header=T,sep="")
####
#geneInfo<-read.delim("last_final_GC_lengths.csv",header=T,sep="")
#a<-a[,c(1,4,3)]
#m<-as.matrix(a[,-1])
#row.names(m)<-a[,1]
#m<-data.frame(m)
#setwd("/Users/asingh2/Desktop/endo_project/counts/")
