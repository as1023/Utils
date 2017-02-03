

https://chitchatr.wordpress.com/2013/06/25/add-error-bars-to-a-plot-in-r/

http://stackoverflow.com/questions/19797846/plot-mean-standard-deviation-standard-error-of-the-mean-and-confidence-interv
http://stackoverflow.com/questions/19797846/plot-mean-standard-deviation-standard-error-of-the-mean-and-confidence-interv

http://www.personality-project.org/r/html/error.bars.html
http://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
http://solvedstack.com/questions/add-error-bars-to-show-standard-deviation-on-a-plot-in-r

https://chitchatr.wordpress.com/2013/06/25/add-error-bars-to-a-plot-in-r/
http://www.cyclismo.org/tutorial/R/confidence.html
http://stackoverflow.com/questions/13224320/r-how-do-i-change-font-style-and-size-on-my-two-x-axes-which-i-have-moved-to-t
http://stackoverflow.com/questions/17551193/r-color-scatter-plot-points-based-on-values
http://jonsullivan.canterburynature.org/?p=332
https://gist.github.com/cdesante/3684833
#####ggplot 
p<-ggplot(data=a,aes(reorder(x=X,-log10(q.val)),y=-log10(q.val)))+ coord_flip() +geom_bar(colour="black", stat="identity",width=.3)+theme_minimal()+ylab(expression(paste('-log10',italic(q),'value')))+ggtitle("Down_GOTerm EvsL")
p <- p + xlab("GO Term")+theme(axis.text = element_text(size = 20),axis.title=element_text(size=18))+theme(axis.ticks = element_line(size = 2))
p<-p+theme(axis.title.x = element_text(face="bold",size=18))
p
dev.off()

