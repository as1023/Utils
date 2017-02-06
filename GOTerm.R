

###load count data in DEseq format and normalize them. convert ENSEMBL name to Entrez id and then to sysmbol
library(gage)
library(gageData)
library(GenomicAlignments)
library(org.Dr.eg.db)
library(pheatmap)

cnts<- counts(DESeq2Table)
sel.rn=rowSums(cnts) != 0
cnts=cnts[sel.rn,]
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
rownames(cnts.norm) <- as.character(mget(rownames(cnts.norm),org.Dr.egENSEMBL2EG,ifnotfound=NA))
rownames(cnts.norm) <- as.character(mget(rownames(cnts.norm), org.Dr.egSYMBOL,ifnotfound=NA))
mymat<-as.matrix(cnts.norm[!grepl('NA', rownames(cnts.norm)), ])

###creat genelist from org.Dr.eg.db

mygo <- as.list(org.Dr.egGO2EG)
mygo <- (mygo[!is.na(mygo)])
t <- mget(names(mygo),GOTERM)
names(mygo) <- as.character(lapply(t,Term))

##### define your data mtarix 
##define control in ref= and define treatment in samp= ... 
gos <- gage(mymat,gsets = mygo,ref=c(1:3),samp=c(4:6),compare ="paired")
gene_set<-sigGeneSet(gos,cutoff=0.01)
out<-(gene_set_1$greater)
write.csv(out,file="go_upregulated.csv")
out<-(gene_set_1$less)
write.csv(out,file="go_downregulated.0001.csv")
#####plot in heatmap 

a<-read.delim("goEvsL_up.0001.csv",sep=",",header=T)
p<-ggplot(data=a[1:15,],aes(reorder(x=X,-log10(q.val)),y=-log10(q.val)))+ coord_flip() +geom_bar(colour="black", stat="identity",width=.3)+theme_minimal()+ylab(expression(paste('-log10',italic(q),'value')))+ggtitle("Up_GOTerm EvsL")
p <- p + xlab("GO Term")+theme(axis.text = element_text(size =13),axis.title=element_text(size=13))+theme(axis.ticks = element_line(size = 2))
p<-p+theme(axis.title.x = element_text(face="bold",size=13))
p
dev.off()





######some webapge 
http://www-huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#gene-ontology-enrichment-analysis
http://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
https://r-forge.r-project.org/scm/viewvc.php/pkg/R/MLPSetUp.R?view=markup&root=mlp&pathrev=21



