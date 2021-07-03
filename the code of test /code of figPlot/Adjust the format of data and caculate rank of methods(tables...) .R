#stable1,2  fig2C

setwd("/Users/xiongkexu/Desktop/双胞项目/测评结果/各数据集结果")
load("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/color.final.robj")
filelist<-Sys.glob("*.csv")
times<-length(filelist)
library(ggplot2)

a<-read.csv(filelist[1],row.names = 1)
auc<-a[,1:1]
pr<-a[,2:2]
for (i in 2:times) {
  auc<-as.data.frame(cbind(auc,read.csv(filelist[i],row.names = 1)[,1:1]))
  pr<-as.data.frame(cbind(pr,read.csv(filelist[i],row.names = 1)[,2:2]))
  a<-a+read.csv(filelist[i],row.names = 1)
}
a<-a/times
a
colnames(auc)<-sub(".csv","",filelist)
colnames(pr)<-sub(".csv","",filelist)
rownames(auc)<-c("bcds","cxds","DoubletFinder","doubletCells","Chord","DoubletDetection","Scrublet","Solo")
rownames(pr)<-c("bcds","cxds","DoubletFinder","doubletCells","Chord","DoubletDetection","Scrublet","Solo")
auc.all<-auc
pr.all<-pr
rm(i)
#调顺序
auc<-auc[order(rownames(auc)),]
pr<-pr[order(rownames(pr)),]

write.csv(round(a,digits = 4),file = "./result/stable1.mean.csv")
write.csv(round(auc,digits = 4),file = "./result/stable2.auroc.csv")
write.csv(round(pr,digits = 4),file = "./result/stable2.pr.csv")

#热图
pdf("fig2c.auc.pdf")
rauc<-apply(auc,2,rank)
pheatmap::pheatmap(auc,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(auc,2,rank)>5,9-apply(auc,2,rank), ""),nrow(auc)),cellwidth = 32,cellheight = 24,angle_col = 45)
dev.off()
pdf("fig2c.pr.pdf")
rpr<-apply(pr,2,rank)
pheatmap::pheatmap(pr,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(pr,2,rank)>5,9-apply(pr,2,rank), ""),nrow(pr)),cellwidth = 32,cellheight = 24,angle_col = 45)
dev.off()






#stable1,2  fig2C

setwd("/Users/xiongkexu/Desktop/双胞项目/测评结果/全方法集成的结果/")
filelist<-Sys.glob("*.csv")
times<-length(filelist)

a<-read.csv(filelist[1],row.names = 1)
auc<-a[,1:1]
pr<-a[,2:2]
for (i in 2:times) {
  auc<-as.data.frame(cbind(auc,read.csv(filelist[i],row.names = 1)[,1:1]))
  pr<-as.data.frame(cbind(pr,read.csv(filelist[i],row.names = 1)[,2:2]))
  a<-a+read.csv(filelist[i],row.names = 1)
}
a<-a/times
a
colnames(auc)<-sub(".csv","",filelist)
colnames(pr)<-sub(".csv","",filelist)
rownames(auc)<-c("bcds","cxds","doubletCells","DoubletFinder","DoubletDetection","Scrublet","Solo","Chord")
rownames(pr)<-c("bcds","cxds","doubletCells","DoubletFinder","DoubletDetection","Scrublet","Solo","Chord")

auc.all<-rbind(auc.all,ChordP=auc["Chord",])
pr.all<-rbind(pr.all,ChordP=pr["Chord",])
rm(i)
#调顺序
auc<-auc[order(rownames(auc)),]
pr<-pr[order(rownames(pr)),]
auc.all<-auc.all[order(rownames(auc.all)),]
pr.all<-pr.all[order(rownames(pr.all)),]

write.csv(round(a,digits = 4),file = "./result/stable1.mean.csv")
write.csv(round(auc,digits = 4),file = "./result/stable2.auroc.csv")
write.csv(round(pr,digits = 4),file = "./result/stable2.pr.csv")

#热图
pdf("fig2c.auc.pdf")
rauc<-apply(auc,2,rank)
pheatmap::pheatmap(auc,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(auc,2,rank)>5,9-apply(auc,2,rank), ""),nrow(auc)),cellwidth = 32,cellheight = 24,angle_col = 45)
dev.off()
pdf("fig2c.pr.pdf")
rpr<-apply(pr,2,rank)
pheatmap::pheatmap(pr,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(pr,2,rank)>5,9-apply(pr,2,rank), ""),nrow(pr)),cellwidth = 32,cellheight = 24,angle_col = 45)
dev.off()


#计算秩排序与整体热图   热图二合一
pdf("./result/auc.all.pdf")
rauc<-apply(auc.all,2,rank)*(-1)+10
pheatmap::pheatmap(auc.all,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(auc.all,2,rank)>6,10-apply(auc.all,2,rank), ""),nrow(auc.all)),cellwidth = 24,cellheight = 18,angle_col = 45)
mat<-as.data.frame(cbind(SDrank=apply(rauc,1,sd),Meanrank=apply(rauc,1,mean),Mean=apply(auc.all,1,mean),Method=rownames(auc.all)))
mat[,1:3]<-apply(mat[,1:3],2,as.numeric)
mat[,1:3]<-round(mat[,1:3],digits = 2)
mat[,4]<-as.factor(mat[,4])
dev.off()

write.csv(mat,file = "./result/auc.all.caculation.csv")
#信息统计
pdf("./result/auc.caculation.pdf",width = 7,height = 3)
p<-ggplot(mat, aes(x=Method, y=SDrank)) + geom_bar(aes(fill=Method),position="dodge", stat="identity")
p+geom_text(aes(label = SDrank), size = 0, hjust = 0.5, vjust = -1, position ="stack")+
  theme_classic()+
#ylab(label = "SDrank")+
#xlab(label = "Method")+
  scale_fill_manual(values = colorlist[[1]])+
  coord_flip()

p<-ggplot(mat, aes(x=Method, y=Meanrank)) + geom_bar(aes(fill=Method),position="dodge", stat="identity")
p+geom_text(aes(label = Meanrank), size = 0, hjust = 0.5, vjust = -1, position ="stack")+
  theme_classic()+
  scale_fill_manual(values = colorlist[[1]])+
  coord_flip()

p<-ggplot(mat, aes(x=Method, y=Mean)) + geom_bar(aes(fill=Method),position="dodge", stat="identity")
p+geom_text(aes(label = Mean), size = 0, hjust = 0.5, vjust = -1, position ="stack")+
  theme_classic()+
  scale_fill_manual(values = colorlist[[1]])+
  coord_flip()
dev.off()

pdf("pr.all.pdf")
rpr<-apply(pr.all,2,rank)*(-1)+10
pheatmap::pheatmap(pr.all,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(pr.all,2,rank)>6,10-apply(pr.all,2,rank), ""),nrow(pr.all)),cellwidth = 32,cellheight = 24,angle_col = 45)
apply(rpr,1,sd)
apply(rpr,1,mean)
apply(pr,1,mean)
dev.off()






#stable1,2  模拟数据集
#setwd("/Users/xiongkexu/Desktop/双胞项目/测评结果/各数据集结果")
setwd("/Users/xiongkexu/Desktop/双胞项目/测评结果/模拟数据集/")
load("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/color.final.robj")
filelist<-Sys.glob("*.csv")
times<-length(filelist)
library(ggplot2)

a<-read.csv(filelist[1],row.names = 1)
auc<-a[,1:1]
pr<-a[,2:2]
for (i in 2:times) {
  auc<-as.data.frame(cbind(auc,read.csv(filelist[i],row.names = 1)[,1:1]))
  pr<-as.data.frame(cbind(pr,read.csv(filelist[i],row.names = 1)[,2:2]))
  a<-a+read.csv(filelist[i],row.names = 1)
}
a<-a/times
a
colnames(auc)<-sub(".csv","",filelist)
colnames(pr)<-sub(".csv","",filelist)
rownames(auc)<-c("bcds","cxds","DoubletFinder","doubletCells","Chord","DoubletDetection","Scrublet","Solo")
rownames(pr)<-c("bcds","cxds","DoubletFinder","doubletCells","Chord","DoubletDetection","Scrublet","Solo")
auc.all<-auc
pr.all<-pr
rm(i)
#调顺序
auc<-auc[order(rownames(auc)),]
pr<-pr[order(rownames(pr)),]

write.csv(round(a,digits = 4),file = "./result/stable1.mean.csv")
write.csv(round(auc,digits = 4),file = "./result/stable2.auroc.csv")
write.csv(round(pr,digits = 4),file = "./result/stable2.pr.csv")

#热图
pdf("fig2c.auc.pdf")
rauc<-apply(auc,2,rank)
pheatmap::pheatmap(auc,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(auc,2,rank)>5,9-apply(auc,2,rank), ""),nrow(auc)),cellwidth = 32,cellheight = 24,angle_col = 45)
dev.off()
pdf("fig2c.pr.pdf")
rpr<-apply(pr,2,rank)
pheatmap::pheatmap(pr,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = matrix(ifelse(apply(pr,2,rank)>5,9-apply(pr,2,rank), ""),nrow(pr)),cellwidth = 32,cellheight = 24,angle_col = 45)
dev.off()



##测试auc95 975 99 
wd1="/Volumes/Transcend/双胞数据/测评/"
wd2=c("A","B","C","batch2.1","batch2.2","HTOdemux8","HTOdemux12")
library("pROC")

pauc80<-data.frame(row.names = wd2,
                   bcds =rep(0,7),
                   cxds =rep(0,7),
                   doubletCells =rep(0,7),
                   DoubletFinder =rep(0,7),
                   DoubletDetection =rep(0,7),
                   Scrublet =rep(0,7),
                   Solo =rep(0,7),
                   ChordP =rep(0,7),
                   Chord =rep(0,7)
)
pauc90<-pauc80
pauc95<-pauc80
pauc975<-pauc80

for (d in 1:length(wd2)) {
  out=wd2[d]
  wd<-paste0(wd1,out,"allmethod")
  mat.all<-read.csv(paste0(wd,"/finalScore.40.csv"))
  mat.old<-read.csv(paste0(sub("allmethod","",wd),"/finalScore.40.csv"))
 
  print(all.equal(mat.all$X,mat.old$X))#检测数据一致性

  setwd("/Users/xiongkexu/Desktop/双胞项目/测评结果/pauc/")
  
  rownames(mat.all)<-mat.all$X
  mat.all<-mat.all[,2:ncol(mat.all)]
  colnames(mat.all)<-c("true","bcds","cxds","doubletCells","DoubletFinder","DoubletDetection","Scrublet","Solo","ChordP")
  mat<-cbind(mat.all,Chord=mat.old$chord)
  
  
  for (i in 2:ncol(mat)) {
    roc<-roc(mat$true,mat[,i],print.auc=TRUE) 
    auc<-auc(roc, # 前面构建的roc1对象
             partial.auc=c(1, 0.025), # 指定计算pAUC的坐标轴范围
             partial.auc.focus="sp", # 根据哪根轴来计算，这里指定特异度范围(x轴)来计算
             partial.auc.correct = F)
    pauc975[out,colnames(mat)[i]]<-as.numeric(auc)
    auc<-auc(roc, # 前面构建的roc1对象
             partial.auc=c(1, 0.05), # 指定计算pAUC的坐标轴范围
             partial.auc.focus="sp", # 根据哪根轴来计算，这里指定特异度范围(x轴)来计算
             partial.auc.correct = F)
    pauc95[out,colnames(mat)[i]]<-as.numeric(auc)
    auc<-auc(roc, # 前面构建的roc1对象
             partial.auc=c(1, 0.1), # 指定计算pAUC的坐标轴范围
             partial.auc.focus="sp", # 根据哪根轴来计算，这里指定特异度范围(x轴)来计算
             partial.auc.correct = F)
    pauc90[out,colnames(mat)[i]]<-as.numeric(auc)
    auc<-auc(roc, # 前面构建的roc1对象
             partial.auc=c(1, 0.2), # 指定计算pAUC的坐标轴范围
             partial.auc.focus="sp", # 根据哪根轴来计算，这里指定特异度范围(x轴)来计算
             partial.auc.correct = F)
    pauc80[out,colnames(mat)[i]]<-as.numeric(auc)
  }
}
write.csv(t(pauc80),file = "pauc80.csv")
write.csv(t(pauc90),file = "pauc90.csv")
write.csv(t(pauc95),file = "pauc95.csv")
write.csv(t(pauc975),file = "pauc975.csv")

pauc<-t(rbind(pauc800=colMeans(pauc80),pauc900=colMeans(pauc90),pauc950=colMeans(pauc95),pauc975=colMeans(pauc975)))

oldauc1<-read.csv("/Users/xiongkexu/Desktop/双胞项目/测评结果/全方法集成的结果/result/stable2.auroc.csv",row.names = 1)
oldauc2<-read.csv("/Users/xiongkexu/Desktop/双胞项目/测评结果/各数据集结果/result/stable2.auroc.csv",row.names = 1)
rownames(oldauc1)[2]<-"ChordP"
oldauc<-rbind(oldauc1["ChordP",],oldauc2)
write.csv(oldauc[colnames(pauc80),],file = "auc.csv")

oldpr1<-read.csv("/Users/xiongkexu/Desktop/双胞项目/测评结果/全方法集成的结果/result/stable2.pr.csv",row.names = 1)
oldpr2<-read.csv("/Users/xiongkexu/Desktop/双胞项目/测评结果/各数据集结果/result/stable2.pr.csv",row.names = 1)
rownames(oldpr1)[2]<-"ChordP"
oldpr<-rbind(oldpr1["ChordP",],oldpr2)
write.csv(oldpr[colnames(pauc80),],file = "pr.csv")

all<-cbind(pauc,auc=rowMeans(oldauc[rownames(pauc),]),pr=rowMeans(oldpr[rownames(pauc),]))
write.csv(all,file = "allaucpr.csv")
