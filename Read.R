
##read file which coloumn is balnk and intriduce NA 
a<-read.delim("test.txt",sep="\t", na.strings=c("", "NA"))
##remove na from column 22
b <- a[!is.na(a[,22]),]
