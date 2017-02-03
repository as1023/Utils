 p<-ggplot(data=a,aes(reorder(x=X,-log10(q.val)),y=-log10(q.val)))+ coord_flip() +geom_bar(colour="black", stat="identity",width=.3)+theme_minimal()+ggtitle("GOTerm H vs L")
 p <- p + xlab("GO Term")+theme(axis.text = element_text(size = 20),axis.title=element_text(size=18,face="bold"))+theme(axis.ticks = element_line(size = 2))
p<-p+theme(axis.title.x = element_text(face="bold",size=18))
