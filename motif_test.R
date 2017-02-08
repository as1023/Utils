library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genome <- BSgenome.Mmusculus.UCSC.mm10
query(MotifDb, "Foxo1")   
pfm.dal80.jaspar <- query(MotifDb,"Foxo1")[[1]]
seqLogo(pfm.dal80.jaspar)

dal1 <- "56458"
up_seqs <- extractUpstreamSeqs(genome, gn, width=1500)
gn <- genes(txdb)[ids.ok]
chromosomal.loc <-transcriptsBy(txdb, by="gene")[dal1]
promoter.dal1 <- getPromoterSeq(chromosomal.loc, genome, upstream=1000, downstream=0)
pcm.dal80.jaspar <- round(100 * pfm.dal80.jaspar)
matchPWM(pcm.dal80.jaspar, unlist(promoter.dal1)[[1]], "90%")
dal80.scertf.hits <- sapply(promoter.dal1, function(pseq) 
                            matchPWM(pcm.dal80.scertf, pseq, min.score="90%"))