
library(DESeq2)
library(biomaRt)
library(Rsamtools)
library(gage)
library(gageData)
library(GenomicAlignments)
library(org.Hs.eg.db)
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
library(edgeR)
library(labdsv)
library('xlsx')
library(statmod)
require(MASS)
library(openxlsx)
library(mclust)
library(Cairo)
CairoFonts(
regular="Arial:style=Regular",
bold="Arial:style=Bold",
italic="Arial:style=Italic",
bolditalic="Arial:style=Bold Italic,BoldItalic",
symbol="Symbol")

####Create a data.frame with the required metadata, i.e., the names of the count files and experimental conditions
#setwd("~/Desktop/LongRNA/human/htseq_out")
 setwd("~/Desktop/LongRNA/human/htseq_out/")
sampleFiles <- c("SN7640334_39846_1_d0_1_STARmapping.htseq_count.txt",
"d0_2_STARmapping.htseq_count.txt","S_3_d0_3_STARmapping.htseq_count.txt",
"d1_2_STARmapping.htseq_count.txt",
"d1_3_STARmapping.htseq_count.txt",
"d1_3_STARmapping.htseq_count.txt",
"d2_1_STARmapping.htseq_count.txt",
"d2_2_STARmapping.htseq_count.txt",
"d2_3_STARmapping.htseq_count.txt",
"d4_1_STARmapping.htseq_count.txt",
"STARmapping.htseq_count.txt",
"d4_3_STARmapping.htseq_count.txt",
"STARmapping.htseq_count.txt",
"d6_2_STARmapping.htseq_count.txt",
"d6_3_STARmapping.htseq_count.txt",
"d8_1_STARmapping.htseq_count.txt",
"d8_2_STARmapping.htseq_count.txt",
"SN7640334_39863_18_d8_3_STARmapping.htseq_count.txt",
"D004270063_39864_19_d14_1_STARmapping.htseq_count.txt",
"D004270063_39865_20_d14_2_STARmapping.htseq_count.txt","D004270063_39866_21_d14_3_STARmapping.htseq_count.txt")
sampleNames <-c("d0_1","d0_2","d0_3","d1_1","d1_2","d1_3","d2_1","d2_2","d2_3","d4_1","d4_2","d4_3","d6_1","d6_2","d6_3","d8_1","d8_2","d8_3","d14_1","d14_2","d14_3")
sampleCondition <- c("d0","d0","d0","d1","d1","d1","d2","d2","d2","d4","d4","d4","d6","d6","d6","d8","d8","d8","d14","d14","d14")
#####Create a CountDataSet object
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = ".",design = ~ condition)
DESeq2Table<-DESeq2Table[rowSums(counts(DESeq2Table)) > 1, ]
#dds<-DESeq(DESeq2Table)
#### Quality control check PCA/samaple cluster, vsd transformation was done to plot PCA. (deseq2 manual -2.2.3)
vsd <- varianceStabilizingTransformation(DESeq2Table, blind=TRUE)
#### PCA plot of the samples with ntop=500 as a defult
CairoPDF(file="Ribo-Samaple_PCA.pdf",width=10/1.54,height=10/1.54)
data<-DESeq2::plotPCA(vsd, intgroup=c("condition"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition))+geom_point(size=3)+ xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+ggtitle("PCA plot Ribo-Samaples")
dev.off()
####Differential regulted gene by Edge R
group<-c("d.0","d.0","d.0","d.1","d.1","d.1","d.2","d.2","d.2","d.4","d.4","d.4","d.6","d.6","d.6","d.8","d.8","d.8","d.14","d.14","d.14")
d<-colnames(assay(DESeq2Table))
y<- DGEList(counts=assay(DESeq2Table), group=group)
keep <- rowSums(cpm(y) > 0.5) >= 2
table(keep)
#keep
y<-y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
options(digits=3)
design <- model.matrix(~ 0+group)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
con <- makeContrasts(d1_vs_d0 = d.1 - d.0, d2_vs_d1 = d.2 - d.1,d4_vs_d2 = d.4-d.2, d6_vs_d4 =d.6 - d.4,d8_vs_d6=d.8-d.6, d14_vs_d8 =d.14 -d.8,levels=design)
res <- glmQLFTest(fit, contrast=con)
######################plot top 50 gene in pheatmap
topDGE <- topTags(res, n = 50)
topDGE <-as.matrix(topDGE$table)
rownames(topDGE) <- as.character(mget(rownames(topDGE),org.Hs.egENSEMBL2EG,ifnotfound=NA))
rownames(topDGE) <- as.character(mget(rownames(topDGE), org.Hs.egSYMBOL,ifnotfound=NA))
CairoPDF(file="Top50 DEG Ribo-genes.pdf",width=12/1.54,height=12/1.54)
pheatmap(topDGE [,1:6],main="Top50 DEG genes(Ribo-",cluster_rows = F,fontsize_row=8)
dev.off()
############
allDGE <- topTags(res , n = nrow(res$table))$table
sigDGE<-subset(allDGE,FDR < 0.01)
####map ensmbel to gene name
sigDGE$ensembl<-rownames(sigDGE)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene","ensembl_transcript_id","refseq_mrna"),filters = "ensembl_gene_id",values = sigDGE$ensembl,mart = ensembl)
idx <- match(sigDGE$ensembl, genemap$ensembl_gene_id)
sigDGE$name <- genemap$hgnc_symbol[idx]
sigDGE$gene <- genemap$entrezgene[idx]
sigDGE$gene_biotype <- genemap$gene_biotype[idx]
#resSig$symbol <- genemap$name_1006[idx]
write.xlsx(data.frame(sigDGE), "DGE_ribo-.xlsx")
