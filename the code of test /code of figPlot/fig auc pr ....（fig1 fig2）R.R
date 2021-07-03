
#setwd("/Volumes/Transcend/双胞数据/测评/A/")
setwd("/Volumes/Transcend/双胞数据/测评/HTOdemux8/")
load("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/color.final.robj")

#fig1 B-------
library(PRROC)
n=1
Fnlist<-c()
dfauc<-data.frame()
dfaupr<-data.frame()
roclist<-list()
prlist<-list()
#数据读入
finalmat<-read.csv("finalScore.40.csv")
#ddc<-read.table("doubletdetection.txt")
#scr<-read.table("scrublet.txt")
#solo<-read.table("solo.txt")

pdf("test2.pdf")
for (i in 3:ncol(finalmat)) {
  group=colnames(finalmat)[i]
  
  true01<-finalmat$true
  true01[which(true01=="Doublet")]<-"1"
  true01[which(true01=="Singlet")]<-"0"
  
  roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =(finalmat[,i]),curve = T)
  if (roc$auc<0.5) {
    finalmat[,i]<--finalmat[,i]
    roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  }
  roclist[[i-2]]<-roc
  plot(roc,main=group)
  pr <- pr.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  prlist[[i-2]]<-pr
  plot(pr,main=group)
  senc<-pr$curve[,1]
  perc<-pr$curve[,2]
  matpr<-cbind(perc,senc)
  
  Fn<-function(n=n,mat=matpr){
    caculateF<-function(x,beta){
      return( ((1+n^2)*(x[1]*x[2]))/(n^2*x[1]+x[2]) )
    }#设置读取每行然后计算Fn的函数，precision，sensitivity作为输入
    
    Fnall<-apply(data.frame(matpr[,1],matpr[,2]), MARGIN=1,FUN =caculateF )#对矩阵
    Fnall[Fnall=="NaN"]<-0  #NA化为0
    return(max(Fnall))
  }
  Fnlist<-c(Fnlist,Fn(n=1,mat=matpr))#记录Fn值
}
dev.off()


color<-as.character(colorlist[[1]][c("bcds","cxds","DoubletFinder","doubletCells","Chord")])

library(Cairo)
CairoPDF("sfig1b.pdf")

for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(roclist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,legend=T)
  }else{
    plot(roclist[[j]],add=TRUE,color=color[j],legend=T)
  }
}

for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(prlist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,)
  }else{
    plot(prlist[[j]],add=TRUE,color=color[j])
  }
}
dev.off()


#fig2 A-------
library(PRROC)
n=1
Fnlist<-c()
dfauc<-data.frame()
dfaupr<-data.frame()
roclist<-list()
prlist<-list()
#数据读入
finalmat<-read.csv("finalScore.40.csv")
ddc<-read.table("doubletdetection.txt")
scr<-read.table("scrublet.txt")
solo<-read.table("solo.txt")
chord_pro<-read.csv("../Aallmethod/finalScore.40.csv")[10]  ####
#chord_pro<-read.csv("../HTOdemux8allmethod//finalScore.40.csv")[10]  
finalmat<-cbind(finalmat,DoubletDetection=ddc,Scrublet=scr,Solo=solo,ChordP=chord_pro)
colnames(finalmat)<-c("X","true","chord","doubletdetection","scrublet","solo","chord_pro")

pdf("test3.pdf")
for (i in 3:ncol(finalmat)) {
  group=colnames(finalmat)[i]
  
  true01<-finalmat$true
  true01[which(true01=="Doublet")]<-"1"
  true01[which(true01=="Singlet")]<-"0"
  
  roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =(finalmat[,i]),curve = T)
  if (roc$auc<0.5) {
    finalmat[,i]<--finalmat[,i]
    roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  }
  roclist[[i-2]]<-roc
  plot(roc,main=group)
  pr <- pr.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  prlist[[i-2]]<-pr
  plot(pr,main=group)
  senc<-pr$curve[,1]
  perc<-pr$curve[,2]
  matpr<-cbind(perc,senc)
  
  Fn<-function(n=n,mat=matpr){
    caculateF<-function(x,beta){
      return( ((1+n^2)*(x[1]*x[2]))/(n^2*x[1]+x[2]) )
    }#设置读取每行然后计算Fn的函数，precision，sensitivity作为输入
    
    Fnall<-apply(data.frame(matpr[,1],matpr[,2]), MARGIN=1,FUN =caculateF )#对矩阵
    Fnall[Fnall=="NaN"]<-0  #NA化为0
    return(max(Fnall))
  }
  Fnlist<-c(Fnlist,Fn(n=1,mat=matpr))#记录Fn值
}
dev.off()

color<-as.character(colorlist[[1]][c("bcds","cxds","DoubletFinder","doubletCells","Chord")])

library(Cairo)
CairoPDF("sfig1B.pdf")

for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(roclist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color3[j],auc.main=FALSE,legend=T)
  }else{
    plot(roclist[[j]],add=TRUE,color=color[j],legend=T)
  }
}

for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(prlist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color3[j],auc.main=FALSE,)
  }else{
    plot(prlist[[j]],add=TRUE,color=color[j])
  }
}
dev.off()





#================================================================================================================

#0316改,fig2全方法的A和HTO8--------

#setwd("/Volumes/Transcend/双胞数据/测评/A/")
setwd("/Volumes/Transcend/双胞数据/测评/HTOdemux8/")
#chord_pro<-read.csv("../Aallmethod/finalScore.40.csv")[10]  ####
chord_pro<-read.csv("../HTOdemux8allmethod//finalScore.40.csv")[10]  
load("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/color.final.robj")

#fig2 A-------
library(PRROC)
n=1
Fnlist<-c()
dfauc<-data.frame()
dfaupr<-data.frame()
roclist<-list()
prlist<-list()
#数据读入
finalmat<-read.csv("finalScore.40.csv")
ddc<-read.table("doubletdetection.txt")
scr<-read.table("scrublet.txt")
solo<-read.table("solo.txt")
finalmat<-cbind(finalmat,DoubletDetection=ddc,Scrublet=scr,Solo=solo,ChordP=chord_pro)
colnames(finalmat)<-c("X","true","bcds","cxds","DoubletFinder","doubletCells","Chord","DoubletDetection","Scrublet","Solo","ChordP")

pdf("test0316.pdf")
for (i in 3:ncol(finalmat)) {
  group=colnames(finalmat)[i]
  
  true01<-finalmat$true
  true01[which(true01=="Doublet")]<-"1"
  true01[which(true01=="Singlet")]<-"0"
  
  roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =(finalmat[,i]),curve = T)
  if (roc$auc<0.5) {
    finalmat[,i]<--finalmat[,i]
    roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  }
  roclist[[i-2]]<-roc
  plot(roc,main=group)
  pr <- pr.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  prlist[[i-2]]<-pr
  plot(pr,main=group)
  senc<-pr$curve[,1]
  perc<-pr$curve[,2]
  matpr<-cbind(perc,senc)
  
  Fn<-function(n=n,mat=matpr){
    caculateF<-function(x,beta){
      return( ((1+n^2)*(x[1]*x[2]))/(n^2*x[1]+x[2]) )
    }#设置读取每行然后计算Fn的函数，precision，sensitivity作为输入
    
    Fnall<-apply(data.frame(matpr[,1],matpr[,2]), MARGIN=1,FUN =caculateF )#对矩阵
    Fnall[Fnall=="NaN"]<-0  #NA化为0
    return(max(Fnall))
  }
  Fnlist<-c(Fnlist,Fn(n=1,mat=matpr))#记录Fn值
}
dev.off()


library(Cairo)
CairoPDF("chord and else0316.pdf")

color<-as.character(colorlist[[1]][colnames(finalmat)[3:11]])
for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(roclist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,legend=T,lwd = 2)
  }else{
    plot(roclist[[j]],add=TRUE,color=color[j],legend=T,lwd = 2)
  }
}

for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(prlist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,)
  }else{
    plot(prlist[[j]],add=TRUE,color=color[j])
  }
}
dev.off()


#fig3 补充两个-------
library(PRROC)
n=1
Fnlist<-c()
dfauc<-data.frame()
dfaupr<-data.frame()
roclist<-list()
prlist<-list()
#数据读入
setwd("/Volumes/Transcend/双胞数据/simulated_pse/1/")
finalmat<-read.csv("finalScore.40.csv")
ddc<-read.table("doubletdetection.txt")
scr<-read.table("scrublet.txt")
solo<-read.table("solo2.txt")
finalmat<-cbind(finalmat,DoubletDetection=ddc,Scrublet=scr,Solo=solo)
colnames(finalmat)<-c("X","true","bcds","cxds","DoubletFinder","doubletCells","Chord","DoubletDetection","Scrublet","Solo")

pdf("test0316.pdf")
for (i in 3:ncol(finalmat)) {
  group=colnames(finalmat)[i]
  
  true01<-finalmat$true
  true01[which(true01=="Doublet")]<-"1"
  true01[which(true01=="Singlet")]<-"0"
  
  roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =(finalmat[,i]),curve = T)
  if (roc$auc<0.5) {
    finalmat[,i]<--finalmat[,i]
    roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  }
  roclist[[i-2]]<-roc
  plot(roc,main=group)
  pr <- pr.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  prlist[[i-2]]<-pr
  plot(pr,main=group)
  senc<-pr$curve[,1]
  perc<-pr$curve[,2]
  matpr<-cbind(perc,senc)
  
  Fn<-function(n=n,mat=matpr){
    caculateF<-function(x,beta){
      return( ((1+n^2)*(x[1]*x[2]))/(n^2*x[1]+x[2]) )
    }#设置读取每行然后计算Fn的函数，precision，sensitivity作为输入
    
    Fnall<-apply(data.frame(matpr[,1],matpr[,2]), MARGIN=1,FUN =caculateF )#对矩阵
    Fnall[Fnall=="NaN"]<-0  #NA化为0
    return(max(Fnall))
  }
  Fnlist<-c(Fnlist,Fn(n=1,mat=matpr))#记录Fn值
}
dev.off()


library(Cairo)
CairoPDF("chord and else0316.pdf")

color<-as.character(colorlist[[1]][colnames(finalmat)[3:10]])
for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(roclist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,legend=T,lwd = 2)
  }else{
    plot(roclist[[j]],add=TRUE,color=color[j],legend=T,lwd = 2)
  }
}

for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(prlist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,)
  }else{
    plot(prlist[[j]],add=TRUE,color=color[j])
  }
}
dev.off()





#数据读入 DEG
setwd("/Volumes/Transcend/双胞数据/simulated data2/")
finalmat<-read.csv("finalScore.40.csv")
ddc<-read.table("doubletdetection.txt")
scr<-read.table("scrublet.txt")
solo<-read.table("solo2.txt")
finalmat<-cbind(finalmat,DoubletDetection=ddc,Scrublet=scr,Solo=solo)
colnames(finalmat)<-c("X","true","bcds","cxds","DoubletFinder","doubletCells","Chord","DoubletDetection","Scrublet","Solo")

pdf("test0316.pdf")
for (i in 3:ncol(finalmat)) {
  group=colnames(finalmat)[i]
  
  true01<-finalmat$true
  true01[which(true01=="Doublet")]<-"1"
  true01[which(true01=="Singlet")]<-"0"
  
  roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =(finalmat[,i]),curve = T)
  if (roc$auc<0.5) {
    finalmat[,i]<--finalmat[,i]
    roc <- roc.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  }
  roclist[[i-2]]<-roc
  plot(roc,main=group)
  pr <- pr.curve( weights.class0 = as.numeric(true01), scores.class0 =finalmat[,i],curve = T)
  prlist[[i-2]]<-pr
  plot(pr,main=group)
  senc<-pr$curve[,1]
  perc<-pr$curve[,2]
  matpr<-cbind(perc,senc)
  
  Fn<-function(n=n,mat=matpr){
    caculateF<-function(x,beta){
      return( ((1+n^2)*(x[1]*x[2]))/(n^2*x[1]+x[2]) )
    }#设置读取每行然后计算Fn的函数，precision，sensitivity作为输入
    
    Fnall<-apply(data.frame(matpr[,1],matpr[,2]), MARGIN=1,FUN =caculateF )#对矩阵
    Fnall[Fnall=="NaN"]<-0  #NA化为0
    return(max(Fnall))
  }
  Fnlist<-c(Fnlist,Fn(n=1,mat=matpr))#记录Fn值
}
dev.off()


library(Cairo)
CairoPDF("chord and else0316.pdf")

color<-as.character(colorlist[[1]][colnames(finalmat)[3:10]])
for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(roclist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,legend=T,lwd = 2)
  }else{
    plot(roclist[[j]],add=TRUE,color=color[j],legend=T,lwd = 2)
  }
}

for (j in 1:(ncol(finalmat)-2)) {
  if (j==1) {
    plot(prlist[[j]],max.plot=TRUE,min.plot=TRUE,rand.plot=TRUE,fill.area=T,color=color[j],auc.main=FALSE,)
  }else{
    plot(prlist[[j]],add=TRUE,color=color[j])
  }
}
dev.off()


