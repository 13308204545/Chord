#setwd("/Volumes/Transcend/双胞数据/测评/HTOdemux8")
setwd("/Volumes/Transcend/双胞数据/测评/A/")
color2 = c("#A6CEE3","#1F78B4","#33A02C","#FFD700","#E31A1C","blue","green","yellow","#FB9A99")
names(color2)=c("bcds","cxds","doubletfinder","doubletcell","chord","doubletdetection","scrublet","solo","chord_nk")

library(Cairo)
library(Seurat)
library(ggplot2)

load("seu2_nooverkill_kmeansclustered.robj")
#提取打分
finalmat<-read.csv("finalScore.40.csv")
ddc<-read.table("doubletdetection.txt")
scr<-read.table("scrublet.txt")
solo<-read.table("solo.txt")
finalmat<-cbind(finalmat,ddc)
finalmat<-cbind(finalmat,scr)
finalmat<-cbind(finalmat,solo)
finalmat$chord<--finalmat$chord
colnames(finalmat)[8:10]<-c("doubletdetection","scrublet","solo")
seu2@meta.data<-cbind(seu2@meta.data[,1:5],finalmat[,2:ncol(finalmat)])
doubletrate<-sum(seu2$true=="Doublet")/length(seu2$true)

seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2, selection.method = "vst", nfeatures = 2000)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2, features = VariableFeatures(object = seu2))
seu2 <- FindNeighbors(seu2, dims = 1:10)
seu2 <- FindClusters(seu2, resolution = 1.2)
seu2 <- RunUMAP(seu2, dims = 1:10)

CairoPDF("sFig1(A).pdf")
alpha <- 0.7
p<-DimPlot(seu2, reduction = "umap",label = F ,pt.size = 0.5,group.by = "true")+
  ggtitle("UMAP") +
  theme_classic(base_size=16) +
  theme(plot.title = element_text(hjust=0.5, size = 22, face='bold'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = 'black'),
        axis.ticks.length=unit(.25, "cm"))

p$layers[[1]]$aes_params$alpha = alpha
p$layers[[1]]$aes_params$stroke = 0
print(p)
print(p+theme(legend.position='n'))
dev.off()

CairoPDF("sFig1A(A).pdf")
alpha <- 0.7
p<-DimPlot(seu2, reduction = "umap",label = F ,pt.size = 0.5,group.by = "true")+
  ggtitle("UMAP") +
  theme_classic(base_size=16) +
  theme(plot.title = element_text(hjust=0.5, size = 22, face='bold'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = 'black'),
        axis.ticks.length=unit(.25, "cm"))

p$layers[[1]]$aes_params$alpha = alpha
p$layers[[1]]$aes_params$stroke = 0
print(p)
print(p+theme(legend.position='n'))
dev.off()

CairoPDF("sFig1B(A).pdf",width = 20,height = 20)
  FeaturePlot(seu2, reduction = "umap",features =colnames(seu2@meta.data)[7:14] ,pt.size = 0.5)
dev.off()




#3.26获取ChordP ==============================
setwd("/Volumes/Transcend/双胞数据/测评/A/")
#setwd("/Volumes/Transcend/双胞数据/测评/HTOdemux8/")
chord_pro<-read.csv("../Aallmethod/finalScore.40.csv")[10]  ####
#chord_pro<-read.csv("../HTOdemux8allmethod//finalScore.40.csv")[10]  
load("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/color.final.robj")

#fig2 A-------
library(PRROC)
library(Seurat)
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
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4//1HTOdemux.R")

##
#HTO8
pbmc.umis <- readRDS("/Volumes/Transcend/双胞数据/测评/HTOdemux8allmethod/pbmc_umi_mtx.rds")
pbmc.htos <- readRDS("/Volumes/Transcend/双胞数据/测评/HTOdemux8allmethod/pbmc_hto_mtx.rds")
seu<-HTO8(pbmc.umis,pbmc.htos,seed=1,out="all",type=8)
x<-table(seu$HTO_classification.global)
doubletrate<-x[1]/(x[1]+x[2])
#A
library(Seurat)
library(Matrix)
counts<-readMM("/Volumes/Transcend/双胞数据/测评//GSE96583_RAW/GSM2560245_A.mat")
barcode<-read.table("/Volumes/Transcend/双胞数据/测评//GSE96583_RAW/GSM2560245_barcodes.tsv")[,1]
genes<-read.table("/Volumes/Transcend/双胞数据/测评//GSE96583_batch1.genes.tsv")[,1]
meta<-read.table("/Volumes/Transcend/双胞数据/测评//GSE96583_batch1.total.tsne.df.tsv",sep = "\t")[barcode,]
rownames(counts)<-genes
colnames(counts)<-barcode
seu<-CreateSeuratObject(counts = counts,)
all.equal(colnames(counts),rownames(meta))
seu$label_scds<-meta$multiplets
seu<-subset(x = seu,subset = label_scds!= "ambs")
doubletrate<-sum(seu$label_scds=="doublet")/length(seu$label_scds)
seu$label_scds[seu$label_scds=="singlet"]<-"Singlet"
seu$label_scds[seu$label_scds=="doublet"]<-"Doublet"
##


seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)
doubletrate<-sum(seu$label_scds=="Doublet")/length(seu$label_scds)

CairoPDF("sFig2(A).pdf")
alpha <- 0.7
p<-DimPlot(seu, reduction = "umap",label = F ,pt.size = 0.5,group.by = "label_scds")+
  ggtitle("UMAP") +
  theme_classic(base_size=16) +
  theme(plot.title = element_text(hjust=0.5, size = 22, face='bold'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = 'black'),
        axis.ticks.length=unit(.25, "cm"))

p$layers[[1]]$aes_params$alpha = alpha
p$layers[[1]]$aes_params$stroke = 0
print(p)
print(p+theme(legend.position='n'))
dev.off()

seu@meta.data<-cbind(seu@meta.data,finalmat[,3:11])#合并

library(Cairo)
CairoPDF("scorePlot.pdf",width = 20,height = 20)
FeaturePlot(seu, reduction = "umap",features =colnames(seu@meta.data)[15:23] ,pt.size = 0.5)
dev.off()

#真阳假阳计算

seu@meta.data$bcds.t<-"Singlet"
seu@meta.data$cxds.t<-"Singlet"
seu@meta.data$doubletCells.t<-"Singlet"
seu@meta.data$DoubletFinder.t<-"Singlet"
seu@meta.data$Chord.t<-"Singlet"       #chord是反的
seu@meta.data$DoubletDetection.t<-"Singlet"
seu@meta.data$Scrublet.t<-"Singlet"
seu@meta.data$Solo.t<-"Singlet"
seu@meta.data$ChordP.t<-"Singlet"      #chord是反的
seu@meta.data[names(sort(x=seu$bcds,decreasing = T))[1:sum(seu$label_scds=="Doublet")],"bcds.t"]<-"Doublet"
seu@meta.data[names(sort(x=seu$cxds,decreasing = T))[1:sum(seu$label_scds=="Doublet")],"cxds.t"]<-"Doublet"
seu@meta.data[names(sort(x=seu$doubletCells,decreasing = T))[1:sum(seu$label_scds=="Doublet")],"doubletCells.t"]<-"Doublet"
seu@meta.data[names(sort(x=seu$DoubletFinder,decreasing = T))[1:sum(seu$label_scds=="Doublet")],"DoubletFinder.t"]<-"Doublet"

seu@meta.data[names(sort(x=seu$Chord,decreasing = F))[1:sum(seu$label_scds=="Doublet")],"Chord.t"]<-"Doublet"

seu@meta.data[names(sort(x=seu$DoubletDetection,decreasing = T))[1:sum(seu$label_scds=="Doublet")],"DoubletDetection.t"]<-"Doublet"
seu@meta.data[names(sort(x=seu$Scrublet,decreasing = T))[1:sum(seu$label_scds=="Doublet")],"Scrublet.t"]<-"Doublet"
seu@meta.data[names(sort(x=seu$Solo,decreasing = T))[1:sum(seu$label_scds=="Doublet")],"Solo.t"]<-"Doublet"

seu@meta.data[names(sort(x=seu$ChordP,decreasing =F))[1:sum(seu$label_scds=="Doublet")],"ChordP.t"]<-"Doublet" 

lable.fpfn<-function(seu=seu,x){
  seu@meta.data[seu$label_scds=="Singlet"&seu@meta.data[,x]=="Doublet",x]<-"FP"
  seu@meta.data[seu$label_scds=="Doublet"&seu@meta.data[,x]=="Singlet",x]<-"FN"
  return(seu)
}
seu<-lable.fpfn(seu,"bcds.t")
seu<-lable.fpfn(seu,"cxds.t")
seu<-lable.fpfn(seu,"doubletCells.t")
seu<-lable.fpfn(seu,"DoubletFinder.t")
seu<-lable.fpfn(seu,"Chord.t")
seu<-lable.fpfn(seu,"DoubletDetection.t")
seu<-lable.fpfn(seu,"Scrublet.t")
seu<-lable.fpfn(seu,"Solo.t")
seu<-lable.fpfn(seu,"ChordP.t")

CairoPDF("resultPlot.pdf",height = 6,width = 6)
alpha <- 0.3
DimPlot(seu,group.by = "bcds.t",cols = c("red","yellow","blue","grey"))
DimPlot(seu,group.by = "bcds.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "cxds.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "doubletCells.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "DoubletFinder.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "Chord.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "DoubletDetection.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "Scrublet.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "Solo.t",cols = c("red","yellow","blue","grey"))+NoLegend()
DimPlot(seu,group.by = "ChordP.t",cols = c("red","yellow","blue","grey"))+NoLegend()  #HTO8的CHORDP 没有FN
#DimPlot(seu2,group.by = "union",cols = c("grey","red"))
#DimPlot(seu2,group.by = "intersection",cols = c("grey","red"))
#p$layers[[1]]$aes_params$alpha = alpha
#p$layers[[1]]$aes_params$stroke = 0
dev.off()


