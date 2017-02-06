
# universe
#resorder is datamatrix 
universe <- rownames(resOrdered)
# load mouse annotation and ID library
biocLite(“org.Mm.eg.db”)
library(org.Mm.eg.db)
# convert gene names to Entrez ID
genemap <- select(org.Mm.eg.db, selected, "ENTREZID", "SYMBOL")
univmap <- select(org.Mm.eg.db, universe, "ENTREZID", "SYMBOL")
# load GO scoring package
biocLite(“GOstats”)
library(GOstats)
# set up analysis
param<- new ("GOHyperGParams", geneIds = genemap, universeGeneIds=univmap, annotation="org.Mm.eg.db",
ontology="BP",pvalueCutoff=0.01, conditional=FALSE, testDirection="over")
# run analysis
hyp<-hyperGTest(param)
# visualize
summary(hyp)
## Select/sort on Pvalue, Count, etc.
###TOP GO
#######Top go analysis
onts = c( "MF", "BP", "CC" )
###take all ensemble id from all the list from deseq2 object 
geneIDs = rownames(norm)
####Select the gene name ..from all gene list to gene interest
geneList <- factor(as.integer (geneIDs %in% d1))
##creat a two level factors
names(geneList) <- geneIDs
tab = as.list(onts)
names(tab) = onts
  for(i in 1:3){
GOdata <- new("topGOdata",ontology=onts[i],allGenes = geneList,nodeSize=5,annot = annFUN.org, mapping="org.Mm.eg.db",ID = "ensembl")
######
resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic <- runTest(GOdata, algorithm = "classic", statistic = "Fisher" )
tab[[i]] <- GenTable(GOdata, Fisher.elim = resultTopGO.elim,Fisher.classic = resultTopGO.classic,orderBy = "Fisher.classic" , topNodes = 25)
}
topGOResults <- rbind.fill(tab)
write.csv(topGOResults, file = "topGOResults.csv")


#######

GOdata <- new("topGOdata",ontology="BP",allGenes = geneList,nodeSize=5,annot = annFUN.org, mapping="org.Mm.eg.db",ID = "ensembl")
######
resultTopGO.elim <- runTest(GOdata, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic <- runTest(GOdata, algorithm = "classic", statistic = "Fisher" )
tab <- GenTable(GOdata, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic,orderBy = "Fisher.classic" , topNodes = 20)


