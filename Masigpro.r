group<-c("d.0","d.0","d.0","d.1","d.1","d.1","d.2","d.2","d.2","d.4","d.4","d.4","d.6","d.6","d.6","d.8","d.8","d.8","d.14","d.14","d.14")
d<-colnames(assay(DESeq2Table))
group<-c("d.0","d.0","d.0","d.1","d.1","d.1","d.2","d.2","d.2","d.4","d.4","d.4","d.6","d.6","d.6","d.8","d.8","d.8","d.14","d.14","d.14")
y<- DGEList(counts=assay(DESeq2Table), group=group)
d<-colnames(assay(DESeq2Table))
a<-read.delim("gene_legth.csv",header=T,sep=",",row.names=1)
y$genes$Length<-a$exon_lengths
norm<-rpkm(y)
#keep <- rowSums(cpm(y) > 0.5) >= 2
#genecount <- cpm(y,normalized.lib.sizes=TRUE)
#design <- model.matrix(~ 0+group)
#colnames(design) <- levels(y$samples$group)
y<- estimateDisp(y, design, robust=TRUE)
#count_data=norm[rowSums(norm > 10) ==ncol(norm), ]
count_data<-norm[rowSums(norm) > 10, ]


#Time <-rep(c(0,1,2,4,6,8,14), each = 3)
#Replicates <- rep(c(1:7), each = 3)
#Group <- rep(1,21)
#ss.edesign <- cbind(Time,Replicates,Group)
#rownames(ss.edesign)<-group

###fetchdesignmetarix file
dd<-read.delim("design.txt",header=T,sep="\t",row.names=1)
d.m <- make.design.matrix(dd,degree=2,time.col=1,repl.col=2,group.cols = c(3:ncol(dd)))
###remove lowcount
#count_data=genecount[rowSums(genecount) > 1, ]
#####use masigprofunction
fit <- p.vector(count_data, d.m, counts=TRUE, theta=1/y$common.dispersion)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
get<-get.siggenes(tstep, vars="all", rsq=0.9)

get$summary
fit$i;
fit$BH.alfa
names(get$sig.genes)
get$sig.profiles
##### extract all gene and pvalue
total_expresion<-cbind(get$sig.genes$sig.profiles,get$sig.genes$sig.pvalues)

#####name symbol and entrez gene
total_expresion$ensembl<-rownames(total_expresion)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene","gene_biotype"),filters = "ensembl_gene_id",values = total_expresion$ensembl,mart = ensembl)


idx <- match(norm$ensembl, genemap$ensembl_gene_id)
total_expresion$name <- genemap$hgnc_symbol[idx]
total_expresion$entrezgene <- genemap$entrezgene[idx]



idx <- match(total_expresion$ensembl, genemap$ensembl_gene_id)
total_expresion$name <- genemap$hgnc_symbol[idx]
total_expresion$entrezgene <- genemap$entrezgene[idx]
#total_expresion$gene_biotype <- genemap$gene_biotype[idx]
write.csv(total_expresion,file="masigpro_sig_gene.csv")
write.xlsx(total_expresion,file="masigpro_sig_gene.xlsx",asTable = TRUE)
