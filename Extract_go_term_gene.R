

b<-as.vector(unlist( mget(ndemo[[""]],org.Mm.egSYMBOL)))
b<-as.data.frame(b)
colnames(b)<-"sym"
t1<-res2[res2$name %in% b$sym,]
t1<-subset(t1,t1$log2FoldChange >0)
k<-data.frame(t1,  row.names=NULL)
t<-k
CairoPDF(file=".pdf",width=15/1.54,height=17/1.54)
p<-ggplot(data=t,aes(reorder(x=t$name,t$log2FoldChange),y=log2FoldChange))+ coord_flip() +geom_bar(colour="black", stat="identity",width=.5)+theme_minimal()+ylab(expression(paste(' ','log2FoldChange')))+ggtitle("")
p <- p + xlab("Gene Name")+theme(axis.text = element_text(size =13),axis.title=element_text(size=13))+theme(axis.ticks = element_line(size = 2))
p
dev.off()


t1<-res1[res1$name %in% b$sym,]
t1<-subset(t1,t1$log2FoldChange >0)
k<-data.frame(t1,  row.names=NULL)
t<-k


CairoPDF(file=".pdf",width=15/1.54,height=17/1.54)
p<-ggplot(data=t,aes(reorder(x=t$name,t$log2FoldChange),y=log2FoldChange))+ coord_flip() +geom_bar(colour="black", stat="identity",width=.5)+theme_minimal()+ylab(expression(paste(' ','log2FoldChange')))+ggtitle("")
p <- p + xlab(" Name")+theme(axis.text = element_text(size =13),axis.title=element_text(size=13))+theme(axis.ticks = element_line(size = 2))
p
dev.off()
