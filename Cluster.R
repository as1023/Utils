##provide a expression matrix
gene_clust<-get$sig.genes$sig.profiles
xyMclust<-Mclust(get$sig.genes$sig.profiles)
gene_clust$CLUST<-xyMclust$classification

####save alll cluster
cluster1 <- names(which(xyMclust$classification==1))
cluster2 <- names(which(xyMclust$classification==2))
cluster3 <- names(which(xyMclust$classification==3))
cluster4 <- names(which(xyMclust$classification==4))
cluster5 <- names(which(xyMclust$classification==5))
cluster6 <- names(which(xyMclust$classification==6))
cluster7 <- names(which(xyMclust$classification==7))
cluster8 <- names(which(xyMclust$classification==8))
cluster9 <- names(which(xyMclust$classification==9))
cluster1<-gene_clust[rownames(gene_clust) %in%cluster1,]
cluster2<-gene_clust[rownames(gene_clust) %in%cluster2,]
cluster3<-gene_clust[rownames(gene_clust) %in%cluster3,]
cluster4<-gene_clust[rownames(gene_clust) %in%cluster4,]
cluster5<-gene_clust[rownames(gene_clust) %in%cluster5,]
cluster6<-gene_clust[rownames(gene_clust) %in%cluster6,]
cluster7<-gene_clust[rownames(gene_clust) %in%cluster7,]
cluster8<-gene_clust[rownames(gene_clust) %in%cluster8,]
cluster9<-gene_clust[rownames(gene_clust) %in%cluster9,]
All<-rbind(cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,cluster7,cluster8,cluster9)

All$ensembl<-rownames(All)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene","gene_biotype"),filters = "ensembl_gene_id",values = All$ensembl,mart = ensembl)
idx <- match(All$ensembl, genemap$ensembl_gene_id)
All$name <- genemap$hgnc_symbol[idx]
All$entrezgene <- genemap$entrezgene[idx]
write.xlsx(All,file="clsuter.xlsx",asTable = TRUE)
###########plot
library(ggplot2)
library(data.table)

CairoPDF(file="All_gene_cluster.pdf",width=15/1.54,height=15/1.54)
mat <- read.delim("All_gene_cluster.csv", sep = ",")
mat$CLUST <- as.character(mat$CLUST)
mat <- mat[,c(2:23,25)]
mat.avg <- t(apply(mat, 1, function(i)
return(c(mean(as.numeric(i[1:3])), mean(as.numeric(i[4:6])), mean(as.numeric(i[7:9])),
mean(as.numeric(i[10:12])), mean(as.numeric(i[13:15])), mean(as.numeric(i[16:18])), mean(as.numeric(i[19:21]))))))
colnames(mat.avg) <- c(0,1,2,4,6,8,14)
ggMat <- melt(data.frame(mat.avg, CLUST = mat$CLUST, name = mat$name))
ggMat$variable <- gsub("X", "", ggMat$variable)
ggMat$variable <- as.numeric(ggMat$variable)
p <- ggplot(ggMat, aes(x=variable, y = value, colour = CLUST, group = name))
p <- p + geom_point(size=1) + geom_line() + facet_wrap( ~ CLUST, nrow = 3, scales = "free_y")
p <- p + stat_summary(aes(group = CLUST), fun.y = mean,   geom = "line", color = "black")+theme_minimal()+ylab("log2Foldchange")+xlab("Time [Day]")
plot(p)
dev.off()


