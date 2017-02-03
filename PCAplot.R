
Skip to content
This repository

    Pull requests
    Issues
    Gist

    @as1023

2
0

    0

dieterich-lab/Leo_zebrafish_downstream_analysis Private
Code
Issues 0
Pull requests 0
Projects 0
Wiki
Pulse
Graphs
Settings
Leo_zebrafish_downstream_analysis/PCA plot
70554b9 on 2 Sep 2016
@as1023 as1023 Rename .gitignore to PCA plot
114 lines (97 sloc) 4.62 KB

library(Cairo)
CairoFonts(
regular="Arial:style=Regular",
bold="Arial:style=Bold",
italic="Arial:style=Italic",
bolditalic="Arial:style=Bold Italic,BoldItalic",
symbol="Symbol")
library(DESeq2)
library(biomaRt)
library(Rsamtools)
library(gage)
library(gageData)
library(GenomicAlignments)
library(org.Dr.eg.db)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(geneplotter)
library(ggplot2)
library(genefilter)
library(annotate)
library(scales)
library(calibrate)
library("Hmisc")
library("GO.db")
library("GOstats")
library(KEGG.db)
library(maSigPro)
library(labdsv)
library('xlsx')
####Create a data.frame with the required metadata, i.e., the names of the count files and experimental conditions
#setwd("/Users/amitsingh/Desktop/leo_project/count_data/zebrafish/")
########load sample from desktop

sampleFiles <- c("SN7640312_35745_0dpi1_STARmapping.htseq_count.txt",
"SN7640312_35746_0dpi2_STARmapping.htseq_count.txt",
"SN7640312_35747_0dpi3_STARmapping.htseq_count.txt","SN7640312_35748_1dpi1_STARmapping.htseq_count.txt",
"SN7640312_35749_1dpi2_STARmapping.htseq_count.txt","SN7640312_35750_1dpi3_STARmapping.htseq_count.txt",
"SN7640312_35751_3dpi1_STARmapping.htseq_count.txt","SN7640312_35752_3dpi2_STARmapping.htseq_count.txt",
"SN7640312_35753_3dpi3_STARmapping.htseq_count.txt",
"SN7640312_35754_7dpi1_STARmapping.htseq_count.txt",
"SN7640312_35756_7dpi3_STARmapping.htseq_count.txt",
"SN7640312_35757_14dpi1_STARmapping.htseq_count.txt",
"SN7640312_35758_14dpi2_STARmapping.htseq_count.txt",
"SN7640312_35759_14dpi3_STARmapping.htseq_count.txt",
"SN7640313_35760_30dpi1_STARmapping.htseq_count.txt",
"SN7640313_35761_30dpi2_STARmapping.htseq_count.txt",
"SN7640313_35762_30dpi3_STARmapping.htseq_count.txt")

sampleNames <- c("d0_1","d0_2","d0_3","d1_1","d1_2","d1_3","d3_1","d3_2","d3_3","d7_1","d7_2","d14_1","d14_2","d14_3","d30_1","d30_2","d30_3")
sampleCondition <- c("d0","d0","d0","d1","d1","d1","d3","d3","d3","d7","d7","d14","d14","d14","d30","d30","d30")
time<-c("0","0","0","1","1","1","3","3","3","7","7","14","14","14","30","30","30")
#####Create a CountDataSet object
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition,time=time)
sampleTable$time <- factor(sampleTable$time, levels=c("0","1","3","7","14","30"))

 
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = ".",design = ~ condition)

#### Quality control check PCA/samaple cluster, vsd transformation was done to plot PCA. (deseq2 manual -2.2.3)
vsd <- varianceStabilizingTransformation(DESeq2Table, blind=TRUE)
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)
#### Plot heatmap to check the correlation between the samples use either of pheatmap of heatmap.2
#pheatmap(mat,cellwidth=10,cellheight=15,cluster_col=T,cex=0.8,cluster_col_method="maximum")

CairoPDF(file="samaple_cluster",width=12/1.54,height=10/1.54)
hmcol <- colorRampPalette(brewer.pal(9, "Reds"))(255)
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

#### PCA plot of the samples with ntop=500 as a defult
CairoPDF(file="Samaple_PCA",width=10/1.54,height=10/1.54)
DESeq2::plotPCA(vsd, intgroup=c("condition"))
dev.off()

#####PCA line trajectories plot in below code 
#cnts<- counts(DESeq2Table)
#sel.rn=rowSums(cnts) != 0
#cnts=cnts[sel.rn,]
#libsizes=colSums(cnts)
#size.factor=libsizes/exp(mean(log(libsizes)))
#cnts.norm=t(t(cnts)/size.factor)
#rownames(cnts.norm) <- as.character(mget(rownames(cnts.norm),org.Dr.egENSEMBL2EG,ifnotfound=NA))
#rownames(cnts.norm) <- as.character(mget(rownames(cnts.norm), org.Dr.egSYMBOL,ifnotfound=NA))
#mymat<-as.matrix(cnts.norm[!grepl('NA', rownames(cnts.norm)), ])
#p5 <- labdsv::pca(t(mymat),dim=4,cor=F)
#mysum <- summary(p5)
#plot(p5$score[,1:2],pch=19,ylab=paste("PC2 (",round(100*mysum[2,2],2),"%)",sep=""),type="n", xlab=paste("PC1 (",round(100*mysum[2,1],2),"%)",sep=""))
#loads <- as.data.frame(p5$loadings[,1:2])
#loads$sym <- as.character(mget(rownames(loads),org.Dr.egSYMBOL))
#loads$len <- sqrt(loads[,1]^2+loads[,2]^2)
#loads <- loads[order(abs(loads[,4]),decreasing=T),]
#num <- 10
#arrows(0,0,200*loads[1:num,1],200*loads[1:num,2],col="grey",length=0.)
#text(200*loads[1:num,1],200*loads[1:num,2],loads[1:num,3],cex=0.5,col="blue")
#text(p5$score[,1],p5$score[,2],rownames(p5$score),pos=3)
#lines(p5$score[1:17,1],p5$score[1:17,2],col="red",type="o",pch=19)
########PCA plot with rlog transformation from DESeq2
#CairoPDF(file="Samaple_PCA",width=10/1.54,height=10/1.54)
#rld <- rlogTransformation(DESeq2Table, blind=T)
#plotPCA(rld, intgroup=c("condition"))
#dev.off()


######Now compute Differential regulated genes 




