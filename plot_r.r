
#### boolean model 
library("BoolNet")
library(colorRamps)
library(pheatmap)
library(labdsv)
library(Cairo)
CairoFonts(
regular="Arial:style=Regular",
bold="Arial:style=Bold",
italic="Arial:style=Italic",
bolditalic="Arial:style=Bold Italic,BoldItalic",
symbol="Symbol"
)



plotHeatmap <- function(n, mat)
{
	refValue = max(c(abs(min(mat)), abs(max(mat))))
	br <- unique(c(seq(-(abs(refValue)), 0, length=6), 0, seq(0,abs(refValue), length=6)))
	pdfH = max(c(4, 0.2*dim(mat)[1]))
	pheatmap(mat, cluster_cols = F, cellwidth = 20, breaks = br, filename = paste(n, ".pdf", sep = ""), height = pdfH, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10))
}


ngf<-loadNetwork("ngf_bool_time2013.txt")
mycol <- matlab.like(12)
mycol[1] <- "#FFFFFF"
#myoder<-c(1:13,58,20,21,16:19,14,15,22:40,59,41:47,48:56,57)
#mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,19),rep(7,8),rep(8,9),rep(11,1))

mycols <- c(rep(1,2),rep(2,4),rep(3,2),rep(4,2),rep(5,14),rep(6,2),rep(7,1))
pheatmap(t(path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8,cluster_col_method="maximum")



### 1st-hour-on 
CairoPDF(file="MAPK_on",width=10/1.54,height=8/1.54)
path1<-getPathToAttractor(ngf,c(1,rep(0,26)))
#out_path1<-path1[,myoder]
pheatmap(t(path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8,cluster_col_method="maximum")
dev.off()

CairoPDF(file="MAPK_inhbition",width=10/1.54,height=8/1.54)
n<-fixGenes(ngf,"MEK",0)
path1<-getPathToAttractor(n,c(1,rep(0,5),0,rep(0,20)))
#ut_path1<-path1[,myoder]
pheatmap(t(path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()


CairoPDF(file="JNK_inhbition",width=10/1.54,height=8/1.54)
n<-fixGenes(ngf,"JNK",0)
path1<-getPathToAttractor(n,c(1,rep(0,5),0,rep(0,20)))
#ut_path1<-path1[,myoder]
pheatmap(t(path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()


CairoPDF(file="PI3k_inhbition",width=10/1.54,height=8/1.54)
n<-fixGenes(ngf,"PI3K",0)
path1<-getPathToAttractor(n,c(1,rep(0,5),0,rep(0,20)))
#ut_path1<-path1[,myoder]
pheatmap(t(path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()







#dev.copy2pdf(file="0-1hr.pdf")

### SERPINE1 on

CairoPDF(file="serpine_on.pdf",width=10/1.54,height=10/1.54)
net1=path1[nrow(path1),] 
net1$SERPINE1<-1
path2<-getPathToAttractor(ngf,as.numeric(net1))
out_path2<-path2[,myoder]
pheatmap(t(out_path2)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
#dev.copy2pdf(file="il6_on_1-3hr.pdf")
dev.off()

### NPY_on 
CairoPDF(file="NPY_on.pdf",width=10/1.54,height=10/1.54)
net2=path2[nrow(path2),]
net2$NPY<-1
net2$SERPINE1<-1
path3<-getPathToAttractor(ngf,as.numeric(net2))
out_path3<-path3[,myoder]
pheatmap(t(out_path3)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()






#######script for inhbition 
CairoPDF(file="MEK_inhbition.pdf",width=10/1.54,height=10/1.54)
#ngf<-loadNetwork("ngf_bool_time2013_inhi.txt")
#mycol <- matlab.like2(12)
#mycol[1] <- "#FFFFFF"
#myoder<-c(1:13,20:22,16:19,14,15,23:38,42,39:41,43:55,56:62)
#mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))
n<-fixGenes(ngf,"MEK",0)
path1<-getPathToAttractor(n,c(1,rep(0,12),0,rep(0,45)))
out_path1<-path1[,myoder]
#dummy1 = c(rep(1,nrow(path1)))
#path1_wi_dummy <- cbind(path1, dummy1)
#out_path1<-path1_wi_dummy[,myoder]
pheatmap(t(out_path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)

dev.off()

####jnk inhbition 

CairoPDF(file="JNK_inhbition.pdf",width=10/1.54,height=10/1.54)
#ngf<-loadNetwork("ngf_bool_time2013_inhi.txt")
#mycol <- matlab.like2(12)
#mycol[1] <- "#FFFFFF"
#myoder<-c(1:13,20:22,16:19,14,15,23:38,42,39:41,43:55,56:62)
#mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))
n<-fixGenes(ngf,"JNK",0)
path1<-getPathToAttractor(n,c(1,rep(0,12),0,rep(0,45)))
out_path1<-path1[,myoder]
pheatmap(t(out_path1)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()


###pi3k inhibtion 


CairoPDF(file="PI3K_inhbition.pdf",width=10/1.54,height=10/1.54)
#ngf<-loadNetwork("ngf_bool_time2013_inhi.txt")
#mycol <- matlab.like2(12)
#mycol[1] <- "#FFFFFF"
#myoder<-c(1:13,20:22,16:19,14,15,23:38,42,39:41,43:55,56:62)
#mycols <- c(rep(1,2),rep(2,11),rep(3,3),rep(4,4),rep(5,2),rep(6,16),rep(7,8),rep(8,9),rep(11,6),rep(12,1))
n<-fixGenes(ngf,"PI3K",0)
path1<-getPathToAttractor(n,c(1,rep(0,12),0,rep(0,45)))
net2=path1[nrow(path1),]
net2$NPY<-1
net2$SERPINE1<-1
path3<-getPathToAttractor(n,as.numeric(net2))
out_path3<-path3[,myoder]
pheatmap(t(out_path3)*mycols,cluster_cols=F,cluster_rows=F,cellwidth=8,color=mycol,cex=0.8)
dev.off()












ngf_uo<-read.csv("ngf+uo126.csv",sep=",",header=T)
ngf<-read.csv("ngf.csv",header=T,sep=",")
ngf_ly<-read.csv("ngf_ly.csv",sep=",",header=TRUE)
ngf_sp<-read.csv("ngf+sp.csv",sep=",",header=T)
#colnames(ngf)<-c("Klf2","Klf4","Btg2","Dusp6","Maff","Jund","Fosl1","Junb","Itga1","Klf5","Klf6","Myc","Pik3ip1","Tieg1","Cited2","Zfp36","Dgka","Prkcd","plat","Npy","Mmp10","Plaur","Mmp3")

colnames(ngf)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")
colnames(ngf_uo)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")
colnames(ngf_sp)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")
colnames(ngf_ly)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")
pdf(file='all.pdf', height=10, width=12)
xaxis=c(0,1,2,4,6,8,12,24,48)

par(mfrow=c(3,3))
limits =  c(-2,14)


for(i in 1:23)
{
	ngf.point=ngf[i,2:10]
	ngf_uo.point=ngf_uo[i,2:10]
	ngf_sp.point=ngf_sp[i,2:10]
	ngf_ly.point=ngf_ly[i,2:10]
	
	plot(xaxis,ngf.point,type="l",ylab="Gene Expression Fold(log2)",xlab="Time (hours)",xlim=c(0,48),main=ngf$gene[i],ylim=limits)
	points(xaxis,ngf.point,pch=18,col="blue",lwd = 5)
	lines(xaxis,ngf.point,col="blue")
	
	
	lines(xaxis,ngf_uo.point,col="red")
	points(xaxis,ngf_uo.point,pch=18,col="red",lwd = 5)
	
	lines(xaxis,ngf_sp.point,col="green")
	points(xaxis,ngf_sp.point,pch=18,col="green",lwd = 5)
	
	lines(xaxis,ngf_ly.point,col="black")
	points(xaxis,ngf_ly.point,pch=18,col="black",lwd = 5)
	
legend("topright",c("NGF","NGF+UO126","NGF+SP","NGF+LY"),fill=c("blue","red","green","black"))
	
}
	
dev.off()
















#simResults <- simulatorDT(CNOlist=cnolist2,model=model,simList=fields4Sim,indices=indexOrig,boolUpdates=10)
#simResults = convert2array(simResults, dim(cnolist2$valueSignals[[1]])[1],length(model$namesSpecies), boolUpdates)
#optimRes <- getFitDT(simResults=simResults,CNOlist=cnolist2,model=model,indexList=indexOrig,boolUpdates=boolUpdates,lowerB=0.5,upperB=10,nInTot=length(which(model$interMat == -1)))
#a<-read.csv("final_data_with_log2fc_4model.csv,sep=",",header=T)
###hil function 

hill<- function(x){
	HillCoef=2
	x^HillCoef/(EC50Data^HillCoef+x^HillCoef)
}

##just tochek 
hill(-0.96900621)
###logistic function

logic <- function(x){
1/(1+exp(-x))
}
###logic(--0.96900621)

#mydata=lapply(list.files(pattern="csv"),function(filename) {
#   dum=read.csv(filename)
#   dum$filename=filename
#   return(dum)
#})

opt<-gaBinaryT1(cnolist2, model, initBstring=init, verbose=TRUE,maxTime=100,relTol=0.01,popSize=100,pMutation=0.5)
pdf(file="sim1hr.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(cnolist2,model,bStrings=list(opt$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()

pdf(file="fit1hr.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt)
dev.off()
#####
pdf(file="fit-model-1hr.pdf",width=18/1.54,height=15/1.54)
plotModel(model, cnolist2, bString=opt$bString)
dev.off()
####plotback model
bs = mapBack(model,pknmodel, opt$bString)
pdf(file="fit-model1_1-hr.pdf",width=18/1.54,height=15/1.54)
plotModel(pknmodel, cnolist2, bs, compressed=model$speciesCompressed)
save(opt,file="res_1h.mat")


o###one time at2hr
t = 2
CNOlistSS =  cnolist2
tIndex = which(CNOlistSS$timeSignals == t)
CNOlistSS$timeSignals = c(0, t)
CNOlistSS$valueSignals = list(t0 =cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex]])
opt1 <- gaBinaryT1(CNOlist = CNOlistSS, model = model,verbose = TRUE,maxTime=100,relTol=0.01,popSize=50,pMutation=0.5)

pdf(file="2hr_sim.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS, model, bStrings = list(opt1$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()



pdf(file="fit2hr.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt1)
dev.off()

#####
pdf(file="fit-model-2hr.pdf",width=18/1.54,height=15/1.54)
plotModel(model,CNOlistSS, bString=opt1$bString)
dev.off()
####plotback model
pdf(file="fit-model2-2hr.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel,opt1$bString)
plotModel(pknmodel, CNOlistSS, bs, compressed=model$speciesCompressed)
dev.off()

####
t = 4
CNOlistSS =  cnolist2
tIndex = which(CNOlistSS$timeSignals == t)
CNOlistSS$timeSignals = c(0, t)
CNOlistSS$valueSignals = list(t0 =cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex]])
opt2 <- gaBinaryT1(CNOlist = CNOlistSS, model = model,verbose = TRUE,maxTime=100,relTol=0.01,popSize=50,pMutation=0.5)
pdf(file="4hr_sim.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS, model, bStrings = list(opt2$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()

pdf(file="4hr_fit.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt2)
dev.off()

#####
pdf(file="fit-model_4hr.pdf",width=18/1.54,height=15/1.54)
plotModel(model,CNOlistSS, bString=opt2$bString,)
dev.off()
####plotback model
pdf(file="fit-model_1_4hr.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel, opt2$bString)
plotModel(pknmodel, CNOlistSS, bs, compressed=model$speciesCompressed)
dev.off()
t = 6
CNOlistSS =  cnolist2
tIndex = which(CNOlistSS$timeSignals == t)
CNOlistSS$timeSignals = c(0, t)
CNOlistSS$valueSignals = list(t0 =cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex]])
opt3 <- gaBinaryT1(CNOlist = CNOlistSS, model = model,verbose = TRUE,maxTime=100,relTol=0.01,popSize=50,pMutation=0.5)

pdf(file="6hr_sim.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS, model, bStrings = list(opt3$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()
pdf(file="6hr_fit.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt3)
dev.off()
#####
pdf(file="6hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model,CNOlistSS, bString=opt3$bString)
dev.off()
####plotback model

pdf(file="6hr_fit_model.1.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel, opt3$bString)
plotModel(pknmodel, CNOlistSS, bs, compressed=model$speciesCompressed)
dev.off()
save(opt3,file="res_6hr.mat")

t = 8
CNOlistSS =  cnolist2
tIndex = which(CNOlistSS$timeSignals == t)
CNOlistSS$timeSignals = c(0, t)
CNOlistSS$valueSignals = list(t0 =cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex]])
opt4 <- gaBinaryT1(CNOlist = CNOlistSS, model = model,verbose = TRUE,maxTime=100,relTol=0.01,popSize=50,pMutation=0.5)
pdf(file="8hr_sim.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS, model, bStrings = list(opt4$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()

pdf(file="8hr_fit.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt4)
dev.off()


#####
pdf(file="8hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model,CNOlistSS, bString=opt4$bString)
dev.off()

####plotback model
pdf(file="8hr_fit_model1.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel, opt4$bString)
plotModel(pknmodel, CNOlistSS, bs, compressed=model$speciesCompressed)
dev.off()
save(opt4,file="res_8hr.mat")

t = 12
CNOlistSS =  cnolist2
tIndex = which(CNOlistSS$timeSignals == t)
CNOlistSS$timeSignals = c(0, t)
CNOlistSS$valueSignals = list(t0 =cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex]])
opt5 <- gaBinaryT1(CNOlist = CNOlistSS, model = model,verbose = FALSE,maxTime=100,relTol=0.01,popSize=50,pMutation=0.5)
pdf(file="12hr_sim.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS, model, bStrings = list(opt5$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()
pdf(file="12hr_fit.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt5)
dev.off()

#####
pdf(file="12hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model,CNOlistSS, bString=opt5$bString)
dev.off()

####plotback model
pdf(file="12hr_fit_model1.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel, opt5$bString)
plotModel(pknmodel, CNOlistSS, bs, compressed=model$speciesCompressed)
dev.off()

save(opt5,file="res_12hr.mat")


#########
t = 24
CNOlistSS =  cnolist2
tIndex = which(CNOlistSS$timeSignals == t)
CNOlistSS$timeSignals = c(0, t)
CNOlistSS$valueSignals = list(t0 =cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex]])
opt6 <- gaBinaryT1(CNOlist = CNOlistSS, model = model,verbose = FALSE,maxTime=100,relTol=0.01,popSize=50,pMutation=0.5)
pdf(file="24hr_fit_model.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS, model, bStrings = list(opt6$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()
pdf(file="24hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt6)
dev.off()
save(opt6,file="res_24hr.mat")
#####
pdf(file="24hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model,CNOlistSS, bString=opt6$bString)
dev.off()
####plotback model

pdf(file="24hr_fit_model1.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel, opt6$bString)
plotModel(pknmodel, CNOlistSS, bs, compressed=model$speciesCompressed)
dev.off()

#######
t = 48
CNOlistSS =  cnolist2
tIndex = which(CNOlistSS$timeSignals == t)
CNOlistSS$timeSignals = c(0, t)
CNOlistSS$valueSignals = list(t0 =cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex]])
opt7 <- gaBinaryT1(CNOlist = CNOlistSS, model = model,verbose = FALSE,maxTime=100,relTol=0.01,popSize=50,pMutation=0.5)
pdf(file="48hr_fit_model.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS, model, bStrings = list(opt7$bString),plotPDF=TRUE,plotParams=list(cex=0.8, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()
pdf(file="48hr_fit.pdf",width=18/1.54,height=15/1.54)
plotFit(optRes=opt7)
dev.off()

#####
pdf(file="48hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model,CNOlistSS, bString=opt7$bString)
dev.off()


####plotback model
pdf(file="48hr_fit_model1.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel, opt7$bString)
plotModel(pknmodel, CNOlistSS, bs, compressed=model$speciesCompressed)
dev.off()
save(opt6,file="res_48hr.mat")


t = c(12, 24)
CNOlistSS2 = cnolist2
tIndex = which(CNOlistSS2$timeSignals == t[1])
tIndex[2] = which(CNOlistSS2$timeSignals == t[2])
# make a new CNOlist with 2 time points
CNOlistSS2$timeSignals = c(0, t)
CNOlistSS2$valueSignals = list(t0 = cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex[1]]], cnolist2$valueSignals[[tIndex[2]]])

#optPart1B1 <- gaBinaryT1(CNOlist=CNOlistSS2, ModelCutCompressExpand,initBstring=initBstring, verbose=TRUE, maxTime=180)


opt2a <- gaBinaryT1(CNOlist = CNOlistSS2, model = model,maxTime =100, verbose = FALSE)
# optimise T2
opt2b <- gaBinaryT2(CNOlist = CNOlistSS2, bStringT1 = opt2a$bString,model = model, maxTime =100, verbose = FALSE)
pdf(file="12-24hr_fit.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS2, model, list(opt2a$bString,opt2b$bString),plotParams=list(cex=0.5, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()

#####
pdf(file="12-24hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model, cnolist2,opt2a$bString,opt2b$bString)
dev.off()
####plotback model
pdf(file="12-24hr_fit_model1.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel,opt2a$bString)
plotModel(pknmodel, cnolist2, bs, compressed=model$speciesCompressed)
dev.off()
save(opt6,file="res_12-24hr.mat")
######steday state with more timepoint

#t = c(2,4,6,8,12,24,48)
#CNOlistSS2 = cnolist2
#tIndex = which(CNOlistSS2$timeSignals == t[1])
#tIndex[2] = which(CNOlistSS2$timeSignals == t[2])
#tIndex[3] = which(CNOlistSS2$timeSignals == t[3])
#tIndex[4] = which(CNOlistSS2$timeSignals == t[4])
#tIndex[5] = which(CNOlistSS2$timeSignals == t[5])
#tIndex[6] = which(CNOlistSS2$timeSignals == t[6])
#tIndex[7] = which(CNOlistSS2$timeSignals == t[7])
#tIndex[8] = which(CNOlistSS2$timeSignals == t[8])

#CNOlistSS2$timeSignals = c(0, t)
#CNOlistSS2$valueSignals = list(t0 = cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex[1]]], cnolist2$valueSignals[[tIndex[2]]],cnolist2$valueSignals[[tIndex[3]]],cnolist2$valueSignals[[tIndex[4]]],cnolist2$valueSignals[[tIndex[5]]],
#cnolist2$valueSignals[[tIndex[6]]],cnolist2$valueSignals[[tIndex[7]]],cnolist2$valueSignals[[tIndex[8]]])
#opt2a <- gaBinaryT1(CNOlist = CNOlistSS2, model = model,maxTime = 60, verbose = FALSE)
#opt2b <- gaBinaryTN(CNOlist = CNOlistSS2, bStrings=list(opt2a$bString),model = model, maxTime = 60, verbose = FALSE)
#cutAndPlot(CNOlistSS2, model,list(opt2a$bString,opt2b$bString),plotParams=list(cex=0.5, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))

####extraction of the opt after 1000 run


#### plot optimized model

#####

t = c(24, 48)
CNOlistSS2 = cnolist2
tIndex = which(CNOlistSS2$timeSignals == t[1])
tIndex[2] = which(CNOlistSS2$timeSignals == t[2])
# make a new CNOlist with 2 time points
CNOlistSS2$timeSignals = c(0, t)
CNOlistSS2$valueSignals = list(t0 = cnolist2$valueSignals[[1]],cnolist2$valueSignals[[tIndex[1]]], cnolist2$valueSignals[[tIndex[2]]])
opt2a <- gaBinaryT1(CNOlist = CNOlistSS2, model = model,maxTime = 100, verbose = FALSE)
# optimise T2
opt2b <- gaBinaryT2(CNOlist = CNOlistSS2, bStringT1 = opt2a$bString,model = model, maxTime = 100, verbose = FALSE)
pdf(file="24-48hr_fit.pdf",width=18/1.54,height=15/1.54)
cutAndPlot(CNOlistSS2, model, list(opt2a$bString,opt2b$bString),plotParams=list(cex=0.5, cmap_scale=0.5, margin=0.2,maxrow=25, margin=0.1, width=20, height=20))
dev.off()

#####
pdf(file="24-48hr_fit_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model, cnolist2,opt2a$bString,opt2b$bString)
dev.off()
####plotback model
pdf(file="24-48hr_fit_model1.pdf",width=18/1.54,height=15/1.54)
bs = mapBack(model,pknmodel,opt2a$bString)
plotModel(pknmodel, cnolist2, bs, compressed=model$speciesCompressed)
dev.off()
save(opt6,file="res_24-48hr.mat")

#tokeep<-c()
#for (i in 1:67)
#{
#	if(b[[i]]=="black")	tokeep<-c(tokeep,names(b)[i])
#}


###plot
#a<-read.csv("time_1hr.csv",header=T,sep=",")
#b<-as.matrix(a[,-1])
#row.names(b)<-a[,1]
#CairoPDF(file="1hr.pdf",width=8/1.54,height=5/1.54)
#pheatmap((b),cluster_cols=F,cluster_rows=F,cellwidth=8,fontsize=8)
#dev.off()

#a<-read.csv("time_2hr.csv",header=T,sep=",")
#b<-as.matrix(a[,-1])
#row.names(b)<-a[,1]
#CairoPDF(file="2hr.pdf",width=8/1.54,height=5/1.54)
#pheatmap((b),cluster_cols=F,cluster_rows=F,cellwidth=8,fontsize=8)
#dev.off()

#a<-read.csv("time_48hr.csv",header=T,sep=",")
#b<-as.matrix(a[,-1])
#row.names(b)<-a[,1]
#CairoPDF(file="48hr.pdf",width=8/1.54,height=5/1.54)
#pheatmap((b),cluster_cols=F,cluster_rows=F,cellwidth=8,fontsize=8)
#dev.off()




CairoPDF(file="ngf.pdf",width=10/1.54,height=10/1.54)
> pheatmap(t(d),cellwidth=8,cellheight=7,cluster_col=F,cex=0.8,cluster_col_method="maximum",scale="row")
> dev.off()





#writeScaffold(modelComprExpanded=model,optimResT1=opt,optimResT2=NA,modelOriginal=pknmodel,CNOlist=cnolist2)
#writeNetwork(modelOriginal=pknmodel,modelComprExpanded=model,optimResT1=opt,optimResT2=NA,CNOlist=cnolist2)
#namesFilesNGF<-list(dataPlot="NGFGraph.pdf",evolFitT1="evolFitToyT1.pdf",evolFitT2=NA,simResultsT1="SimResultsT1_1.pdf",
#simResultsT2=NA,scaffold="Scaffold.sif",scaffoldDot="Scaffold.dot",tscaffold="TimesScaffold.EA", wscaffold="weightsScaffold.EA", PKN="PKN.sif",PKNdot="PKN.dot",wPKN="TimesPKN.EA", nPKN="nodesPKN.NA")
#writeReport(modelOriginal=pknmodel,modelOpt=model,optimResT1=opt,optimResT2=NA,CNOlist=cnolist2,directory="bestfit",namesFiles=namesFilesNGF,namesData=list(CNOlist="cnolist2",model="pknmodel"))


######ODEmodeling Rod..trainng model with all data 

model <- preprocessing(cnolist2, pknmodel,expansion=TRUE,compression=FALSE,verbose=FALSE,cutNONC=TRUE)
res<-residualError(cnolist2)
init<-rep(1,length(model$reacID))
opt1 <- gaBinaryDT(CNOlist=cnolist2, model=model, initBstring=init,verbose=FALSE, boolUpdates=10, maxTime=1000, lowerB=0.2, upperB=10,maxSizeHashTable=10000,popSize=50, maxGens=1000,relTol=0.01)
pdf(file="all_timepoint_model.pdf",width=25/1.54,height=15/1.54)
cutAndPlotResultsDT(model=model,CNOlist=cnolist2, bString=opt1$bString,plotPDF=FALSE,boolUpdates=10,lowerB=0.2,upperB=10)
dev.off()
###plotback to get optimized model
pdf(file="all_timepoint_opt_model.pdf",width=18/1.54,height=15/1.54)
plotModel(model, cnolist2, bString=opt1$bString)
dev.off()

###rode model fitting

#pknmodel=readSIF("PKN-netwrok.sif")
#a<-load("cnolist.mat")



#initial_pars=createLBodeContPars(pknmodel, LB_n = 1, LB_k = 0.1,
#LB_tau = 0.001, UB_n = 5, UB_k = 0.6, UB_tau = 10, random = TRUE)
#Visualize initial solution
#simulatedData=plotLBodeFitness(cnolist2, pknmodel,initial_pars)
#paramsGA = defaultParametersGA()
#paramsGA$maxStepSize = 1
#paramsGA$popSize = 50
#paramsGA$iter = 100
#paramsGA$transfer_function = 3
#opt_pars=parEstimationLBode(cnolist2,pknmodel,ode_parameters=initial_pars,paramsGA=paramsGA)
#Visualize fitted solution
#simulatedData=plotLBodeFitness(cnolist, pknmodel,ode_parameters=opt_pars)


####
#nmydata <- lapply(mydata, function(x) {c(x, rep(NA, max.len - length(x)))})
#mat <- do.call(rbind, corrected.list)
###
d<-read.csv("ngf.csv",header=T,sep=",")
d<-t(d)
colnames(ngf)<-c("1h","2h","4h","6h","8h","12h","24h","48h")
d<-t(d)
phe
colnames(t)<-c("Klf2","Klf4","Btg2","Dusp6","Maff","Jund","Fosl1","Junb","Itga1","Klf5","Klf6","Myc","Tieg1","Cited2","Zfp36","Dgka","Prkcd","Npy","Mmp10")

colnames(ngf)<-c("Klf2","Klf4","Btg2","Dusp6","Maff","Jund","Fosl1","Junb","Itga1","Klf5","Klf6","Myc","Pik3ip1","Tieg1","Cited2","Zfp36","Dgka","Prkcd","Npy","Mmp10")


###
pdf("uo126+ngf.pdf",height=20,width=7)


#######

ngf_uo<-read.csv("ngf_uo.csv",sep=",",header=T)
ngf<-read.csv("ngf.csv",header=T,sep=",")
ngf_ly<-read.csv("ngf_ly.csv",sep=",",header=T)
ngf_sp<-read.csv("ngf_sp.csv",sep=",",header=T)

colnames(ngf)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")
colnames(ngf_uo)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")
colnames(ngf_sp)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")
colnames(ngf_ly)<-c("gene","0h","1h","2h","4h","6h","8h","12h","24h","48h")

pdf(file="peturbed_1_1mek_1_1.pdf", height=10, width=12)
xaxis=c(0,1,2,4,6,8,12,24,48)
#ticks = c(0,1,2,4,6,8,12,24,48);

par(mfrow=c(3,3))

for(i in 1:23){
ngf.point=ngf[i,2:10]
ngf_uo.point=ngf_uo[i,2:10]
ngf_sp.point=ngf_sp[i,2:10]
ngf_ly.point=ngf_ly[i,2:10]

plot(xaxis,ngf.point,type="l",main=ngf$gene[i],ylim=c(-2,13))
points(xaxis,ngf.point,pch=18,col="blue",cex=2)
lines(xaxis,ngf.point,col="blue",lwd =3)
#axis(side = 1, at = ticks);


lines(xaxis,ngf_uo.point,col="red",type="l",lwd =3)
points(xaxis,ngf_uo.point,pch=18,cex=2,col="red")

lines(xaxis,ngf_sp.point,col="black",type="l",lwd =3)
points(xaxis,ngf_sp.point,pch=18,cex=2,col="black")

#lines(xaxis,ngf_ly.point,col="green",type="l",lwd =3)
#points(xaxis,ngf_ly.point,pch=18,cex=2,col="green")

legend("topleft",c("NGF","NGF+UO126","NGF+SP"),fill=c("blue","red","black"))	

}

dev.off()


