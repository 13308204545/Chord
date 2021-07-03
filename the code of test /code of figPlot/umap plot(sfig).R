library(Seurat)
library(ggplot2)
setwd("/Users/xiongkexu/Desktop/双胞项目/论文写作/fig1/")
#
####data demuxlet A----
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

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)
#倒入打分值
finalmat<-read.csv("/Volumes/Transcend/双胞数据/测评/A/finalScore.40.csv")
seu$Chord<-finalmat$chord
##print
pdf("A.pdf")
FeaturePlot(seu,features = "Chord",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)+NoLegend()
FeaturePlot(seu,features = "Chord",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)+NoLegend()
dev.off()




####data HTO8----
library(Seurat)
library(Matrix)
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4//1HTOdemux.R")

pbmc.umis <- readRDS("/Volumes/Transcend/双胞数据/测评/HTOdemux8allmethod/pbmc_umi_mtx.rds")
pbmc.htos <- readRDS("/Volumes/Transcend/双胞数据/测评/HTOdemux8allmethod/pbmc_hto_mtx.rds")
seu<-HTO8(pbmc.umis,pbmc.htos,seed=1,out="all",type=8)
x<-table(seu$HTO_classification.global)
doubletrate<-x[1]/(x[1]+x[2])

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)

#倒入打分值
finalmat<-read.csv("/Volumes/Transcend/双胞数据/测评/HTOdemux8/finalScore.40.csv")
seu$Chord<-finalmat$chord
##
pdf("HTO8.pdf")
FeaturePlot(seu,features = "Chord",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)+NoLegend()
FeaturePlot(seu,features = "Chord",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)+NoLegend()
dev.off()


####data DEG----
library(Seurat)
library(Matrix)

load("/Volumes/Transcend/双胞数据/simulated data2/seu.robj")

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)

#倒入打分值
finalmat<-read.csv("/Volumes/Transcend/双胞数据/simulated data2/finalScore.40.csv")
seu$Chord<-finalmat$chord
pdf("DEG.pdf")
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)+NoLegend()
dev.off()

load("/Volumes/Transcend/双胞数据/simulated_pse/1/seu.robj")

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)
pdf("pse.pdf")
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)
DimPlot(seu,group.by = "label_scds",cols =c("#EE3B3B","lightgray"),pt.size = 0.1)+NoLegend()
dev.off()




