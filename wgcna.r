####wgcna

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
library(org.Mm.eg.db)
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
library('xlsx')
library(doParallel)
library(foreach)
library(dynamicTreeCut)
library(cluster)
library(flashClust)
library(Hmisc)
library(reshape)

####Load the count data
sampleFiles<- c("E1.htseq_count.txt",
"E2.htseq_count.txt",
"E3.htseq_count.txt",
"L1.htseq_count.txt",
"L2.htseq_count.txt",
"L3.htseq_count.txt",
"H1.htseq_count.txt",
"H2.htseq_count.txt",
"H3.htseq_count.txt")
sampleNames <- c("E1","E2","E3","L1","L2","L3","H1","H2","H3")
sampleCondition <- c("E","E","E","L","L","L","H","H","H")
####following command to Make Deseq2 readable format
#sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
#DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = ".",design = ~ condition)
#cnts<-GeneCounts
#sel.rn=rowSums(cnts) != 0
#cnts=cnts[sel.rn,]
#dim(cnts)
#libsizes=colSums(cnts)
#size.factor=libsizes/exp(mean(log(libsizes)))
#cnts.norm=t(t(cnts)/size.factor)
#range(cnts.norm)
#cnts.norm=log2(cnts.norm+8)
#GeneCounts[filter,]))
#filter <- apply(GeneCounts,1,function(x) mean(x)>10)


######
#cnts.norm=log2(cnts.norm+1)
#dim(cnts.norm)

vsd <- varianceStabilizingTransformation(DESeq2Table, blind=TRUE)
Pvars <- rowVars(assay(vsd))
ntop=6500
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,length(Pvars)))]


trans.oed<-(assay(vsd)[select, ])
head(trans.oed)
gene.names=rownames(trans.oed)
trans.oed=t(trans.oed)
datExpr=trans.oed
dim(datExpr)
SubGeneNames=gene.names

#gene.names=rownames(cnts.norm)
#trans.oed=t(cnts.norm)
#datExpr=trans.oed
#dim(datExpr)
#SubGeneNames=gene.names
powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
#Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower = 12;
#calclute the adjacency matrix
adj= adjacency(datExpr,type = "unsigned", power = softPower);
#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower);
colnames(TOM) =rownames(TOM) =SubGeneNames
dissTOM=1-TOM

##Module detection
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);
minModuleSize = 10;
# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
####
CairoPDF(file="cluste.pdf",width=12/1.54,height=10/1.54)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)
colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")



