library(labdsv)
library(Cairo)
library(org.Mm.eg.db)
CairoFonts(
regular="Arial:style=Regular",
bold="Arial:style=Bold",
italic="Arial:style=Italic",
bolditalic="Arial:style=Bold Italic,BoldItalic",
symbol="Symbol"
)

###trajectory plot
#####combine all your normalized data
mymat <- cbind(t[,4:6],t[,1:3])
##name esemble id to gene symbol
mysym <- mget(rownames(mymat),org.Mm.egSYMBOL,ifnotfound=NA)
###compute pca
p5 <- labdsv::pca(t(mymat),dim=4,cor=F)
mysum <- summary(p5)
##scince I was ploting two type data thats why I uesd two color but for your case would be one color
plot(p5$score[,1:2],pch=19,col=c(rep("green",3),rep("red",3)),
#xlim=c(-25,25),ylim=c(-25,25),type="n",
ylab=paste("PC2 (",round(100*mysum[2,2],2),"%)",sep=""),xlab=paste("PC1 (",round(100*mysum[2,1],2),"%)",sep=""))
##loading annotation thats why you neeed again annotation
mymat <- cbind(t[,4:6],t[,1:3])
mysym <- mget(rownames(mymat),org.Mm.egSYMBOL,ifnotfound=NA)
p5 <- labdsv::pca(t(mymat),dim=4,cor=F)
mysum <- summary(p5)
plot(p5$score[,1:2],pch=19,col=c(rep("green",3),rep("red",3)),ylab=paste("PC2 (",round(100*mysum[2,2],2),"%)",sep=""),xlab=paste("PC1 (",round(100*mysum[2,1],2),"%)",sep=""))
loads <- as.data.frame(p5$loadings[,1:2])
loads$sym <- as.character(mget(rownames(loads),org.Mm.egSYMBOL,ifnotfound=NA))
loads$len <- sqrt(loads[,1]^2+loads[,2]^2)
loads <- loads[order(abs(loads[,4]),decreasing=T),]
num <- 10
arrows(0,0,200*loads[1:num,1],200*loads[1:num,2],col="grey",length=0.)
text(200*loads[1:num,1],200*loads[1:num,2],loads[1:num,3],cex=0.5,col="blue")
text(p5$score[,1],p5$score[,2],rownames(p5$score),pos=3)
##the below code is for specifying the line for contol and treatment
lines(p5$score[1:3,1],p5$score[1:3,2],col="red",type="o",pch=19)
lines(p5$score[4:6,1],p5$score[4:6,2],col="green",type="o",pch=19)
legend("topleft",c("E","L"),fill=c("red","green"),inset=0.02)