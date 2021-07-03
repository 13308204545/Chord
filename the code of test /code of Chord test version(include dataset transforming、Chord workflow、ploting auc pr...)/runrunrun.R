##V4.0 
#添加两种新的方法
#添加组合细胞比例随机机制（高斯分布）
#优化出图

###########======仿真实数据流程------
#runrunrun!
setwd("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/")
source("./1HTOdemux.R")
source("./2SingleCellExperiment_INPUT.R")
source("./3cxds boosting.R")
source("./3DBLFinder.R")
source("./4testauc.R")
source("./5adboosting.R")
source("./1overkillDB.R")
source("./3scran.R")

####data demuxlet A----
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/Aallmethod/")
counts<-readMM("../GSE96583_RAW/GSM2560245_A.mat")
barcode<-read.table("../GSE96583_RAW/GSM2560245_barcodes.tsv")[,1]
genes<-read.table("../GSE96583_batch1.genes.tsv")[,1]
meta<-read.table("../GSE96583_batch1.total.tsne.df.tsv",sep = "\t")[barcode,]
rownames(counts)<-genes
colnames(counts)<-barcode
seu<-CreateSeuratObject(counts = counts,)
all.equal(colnames(counts),rownames(meta))
seu$label_scds<-meta$multiplets
seu<-subset(x = seu,subset = label_scds!= "ambs")
doubletrate<-sum(seu$label_scds=="doublet")/length(seu$label_scds)
seu$label_scds[seu$label_scds=="singlet"]<-"Singlet"
seu$label_scds[seu$label_scds=="doublet"]<-"Doublet"

####data demuxlet B----
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/Ballmethod/")
counts<-readMM("../GSE96583_RAW/GSM2560246_B.mat")
barcode<-read.table("../GSE96583_RAW/GSM2560246_barcodes.tsv")[,1]
genes<-read.table("../GSE96583_batch1.genes.tsv")[,1]
meta<-read.table("../GSE96583_batch1.total.tsne.df.tsv",sep = "\t")[barcode,]
rownames(counts)<-genes
colnames(counts)<-barcode
seu<-CreateSeuratObject(counts = counts,)
all.equal(colnames(counts),rownames(meta))
seu$label_scds<-meta$multiplets
seu<-subset(x = seu,subset = label_scds!= "ambs")
doubletrate<-sum(seu$label_scds=="doublet")/length(seu$label_scds)
seu$label_scds[seu$label_scds=="singlet"]<-"Singlet"
seu$label_scds[seu$label_scds=="doublet"]<-"Doublet"

####data demuxlet C
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/Callmethod/")
counts<-readMM("../GSE96583_RAW/GSM2560247_C.mat")
barcode<-read.table("../GSE96583_RAW/GSM2560247_barcodes.tsv")[,1]
genes<-read.table("../GSE96583_batch1.genes.tsv")[,1]
meta<-read.table("../GSE96583_batch1.total.tsne.df.tsv",sep = "\t")[barcode,]
rownames(counts)<-genes
colnames(counts)<-barcode
seu<-CreateSeuratObject(counts = counts,)
all.equal(colnames(counts),rownames(meta))
seu$label_scds<-meta$multiplets
seu<-subset(x = seu,subset = label_scds!= "ambs")
doubletrate<-sum(seu$label_scds=="doublet")/length(seu$label_scds)
seu$label_scds[seu$label_scds=="singlet"]<-"Singlet"
seu$label_scds[seu$label_scds=="doublet"]<-"Doublet"
####data demuxlet 2.1bitch2------
#
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/batch2.1allmethod/")
counts<-as.matrix(readMM("../GSE96583_RAW/GSM2560248_2.1.mtx"))
barcode<-read.table("../GSE96583_RAW/GSM2560248_barcodes.tsv")[,1]
genes<-read.table("../GSE96583_batch2.genes.tsv")[,1]
meta<-read.table("../GSE96583_batch2.total.tsne.df.tsv",sep = "\t")[barcode,]

rownames(counts)<-genes
colnames(counts)<-barcode
seu<-CreateSeuratObject(counts = counts,)
all.equal(colnames(counts),rownames(meta))
seu$label_scds<-meta$multiplets
seu<-subset(x = seu,subset = label_scds!= "ambs")
doubletrate<-sum(seu$label_scds=="doublet")/length(seu$label_scds)
seu$label_scds[seu$label_scds=="singlet"]<-"Singlet"
seu$label_scds[seu$label_scds=="doublet"]<-"Doublet"

####data demuxlet 2.1bitch3------
#
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/batch2.1allmethod_allkill/")
counts<-as.matrix(readMM("../GSE96583_RAW/GSM2560248_2.1.mtx"))
barcode<-read.table("../GSE96583_RAW/GSM2560248_barcodes.tsv")[,1]
genes<-read.table("../GSE96583_batch2.genes.tsv")[,1]
meta<-read.table("../GSE96583_batch2.total.tsne.df.tsv",sep = "\t")[barcode,]

rownames(counts)<-genes
colnames(counts)<-barcode
seu<-CreateSeuratObject(counts = counts,)
all.equal(colnames(counts),rownames(meta))
seu$label_scds<-meta$multiplets
seu<-subset(x = seu,subset = label_scds!= "ambs")
doubletrate<-sum(seu$label_scds=="doublet")/length(seu$label_scds)
seu$label_scds[seu$label_scds=="singlet"]<-"Singlet"
seu$label_scds[seu$label_scds=="doublet"]<-"Doublet"
####data demuxlet 2.2bitch------
#
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/batch2.2allmethod/")
counts<-as.matrix(readMM("../GSE96583_RAW/GSM2560249_2.2.mtx"))
barcode<-read.table("../GSE96583_RAW/GSM2560249_barcodes.tsv")[,1]
genes<-read.table("../GSE96583_batch2.genes.tsv")[,1]
meta<-read.table("../GSE96583_batch2.total.tsne.df.tsv",sep = "\t")[barcode,]

rownames(counts)<-genes
colnames(counts)<-barcode
seu<-CreateSeuratObject(counts = counts,)
all.equal(colnames(counts),rownames(meta))
seu$label_scds<-meta$multiplets
seu<-subset(x = seu,subset = label_scds!= "ambs")
doubletrate<-sum(seu$label_scds=="doublet")/length(seu$label_scds)
seu$label_scds[seu$label_scds=="singlet"]<-"Singlet"
seu$label_scds[seu$label_scds=="doublet"]<-"Doublet"
####data HTO8 all------
#设置数据集的预期双胞率222222
setwd("/Volumes/Transcend/双胞数据/测评/HTOdemux8allmethod")
pbmc.umis <- readRDS("./pbmc_umi_mtx.rds")
pbmc.htos <- readRDS("./pbmc_hto_mtx.rds")
seu<-HTO8(pbmc.umis,pbmc.htos,seed=1,out="all",type=8)
x<-table(seu$HTO_classification.global)
doubletrate<-x[1]/(x[1]+x[2])

####data HTO8 all kill------
#设置数据集的预期双胞率222222
setwd("/Volumes/Transcend/双胞数据/测评/HTOdemux8allmethodkill/")
pbmc.umis <- readRDS("./pbmc_umi_mtx.rds")
pbmc.htos <- readRDS("./pbmc_hto_mtx.rds")
seu<-HTO8(pbmc.umis,pbmc.htos,seed=1,out="all",type=8)
x<-table(seu$HTO_classification.global)
doubletrate<-x[1]/(x[1]+x[2])

#1data HTO12 all
setwd("/Volumes/Transcend/双胞数据/测评/HTOdemux12allmethodkill/")
pbmc.umis <- readRDS("./hto12_umi_mtx.rds")
pbmc.htos <- readRDS("./hto12_hto_mtx.rds")
seu<-HTO8(pbmc.umis,pbmc.htos,seed=1,out="all",type=12)
x<-table(seu$HTO_classification.global)
doubletrate<-x[1]/(x[1]+x[2])

####data SOLOhm2 18940gene 21179cells------ 
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/GSE140262/")
meta<-read.csv("GSE140262_kidney2metadata.csv")
rownames(meta)<-meta$index
gene<-read.csv("GSE140262_kidney2_gene_data.csv")
counts<-readMM("./GSE140262_RAW/GSM4158565_kidney2_RNA.mtx")

rownames(counts)<-meta$index
colnames(counts)<-gene$index

seu<-CreateSeuratObject(counts = t(counts))
seu@meta.data<-meta
seu$label_scds<-"Singlet"
seu$label_scds[seu$Category_with_clustering%in%"Doublet"]<-"Doublet"
rm(counts)
doubletrate=sum(seu$label_scds=="Doublet")/ncol(seu)

####data SOLOhm2 18940gene 21179cells------ 
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/GSE140262allmethod/")
meta<-read.csv("GSE140262_kidney2metadata.csv")
rownames(meta)<-meta$index
gene<-read.csv("GSE140262_kidney2_gene_data.csv")
counts<-readMM("./GSE140262_RAW/GSM4158565_kidney2_RNA.mtx")

rownames(counts)<-meta$index
colnames(counts)<-gene$index

seu<-CreateSeuratObject(counts = t(counts))
seu@meta.data<-meta
seu$label_scds<-"Singlet"
seu$label_scds[seu$Category_with_clustering%in%"Doublet"]<-"Doublet"
rm(counts)
doubletrate=sum(seu$label_scds=="Doublet")/ncol(seu)

####data SOLOhm1 18940gene 21179cells------ 
library(Seurat)
library(Matrix)
setwd("/Volumes/Transcend/双胞数据/测评/GSE140262k1allmethod")
meta<-read.csv("kidney1metadata.csv")
rownames(meta)<-meta$index
gene<-read.csv("GSE140262_kidney1_gene_data.csv")
counts<-readMM("./GSE140262_RAW/GSM4158563_kidney1_RNA.mtx")

rownames(counts)<-meta$index
colnames(counts)<-gene$index

seu<-CreateSeuratObject(counts = t(counts))
seu@meta.data<-meta
seu$label_scds<-"Singlet"
seu$label_scds[seu$Category_with_clustering%in%"Doublet"]<-"Doublet"
rm(counts)
doubletrate=sum(seu$label_scds=="Doublet")/ncol(seu)

##############run----

mfinal=40
k=20
overkill=T
overkillrate=1


library(reticulate)
library(scater)
library(loomR)
library(Seurat)
library(patchwork)

sce<-creatSCE(seu=seu)

sce<-scds(sce=sce)

seu<-scranDB(seu=seu)

seu<-DBF(seu=seu,ground_truth = F,doubletrate=doubletrate)

mattrain<-testroc(seu=seu,sce=sce,outname = "train")#第一次测试,只对有实验标签的使用
write.csv(mattrain,file="realscore.csv")

#seu2<-overkillDB2(seu=seu,sce=sce,doubletrate=doubletrate,seed=1,out="all",k=k,overkill=overkill,overkillrate=overkillrate,add=add)
#save(seu2,file = "1.22alloverkill.robj")
seu2<-overkillDB2(seu=seu,sce=sce,doubletrate=doubletrate,seed=1,out="all",k=k,overkill=overkill,overkillrate=overkillrate)#V3

doubletrate2=sum(seu2$label_scds=="Doublet")/ncol(seu2)

sce2<-creatSCE(seu=seu2)

sce2<-scds(sce=sce2)

seu2<-scranDB(seu=seu2)

seu2<-DBF(seu=seu2,ground_truth = F,doubletrate=doubletrate2)

save(seu2,file="seu_overkilled_run.robj")
save(sce2,file="sce_overkilled_run.robj")
load("seu_overkilled_run.robj")
load("sce_overkilled_run.robj")


#准备处理另外三种
mx<-seu2@assays$RNA@counts
write.csv(mx,file = "counts2.csv")
use_condaenv("base")
Sys.setenv(RETICULATE_PYTHON = '/anaconda3/bin/python')
use_python('/anaconda3/bin/python',required = T)
source_python("/Users/xiongkexu/Desktop/python test/1x.py")
source_python("/Users/xiongkexu/Desktop/python test/2x.py")
as.loom(x=seu2, filename = "seu22.loom", verbose = FALSE,overwrite = T)

source_python("/Users/xiongkexu/Desktop/python test/3x.py")

#
load("seu_overkilled_run.robj")
load("sce_overkilled_run.robj")
#load("seu3_old and new DB.robj")
#seu2<-seu3[,seu3$overkill=="Singlet"]#改动过
#rm(seu3)#改动过1.20

#处理完solo之后，整合所有要导入的数据为add2
source_python("/Users/xiongkexu/Desktop/python test/4(treat npy)x.py")
solo<-read.table("solo2.txt")
ddc<-read.table("doubletdetection2.txt")
scr<-read.table("scrublet2.txt")
add2<-cbind(ddc,scr,solo)
colnames(add2)<-c("doubletdetection","scrublet","solo")
rownames(add2)<-colnames(seu2)
#
mattrain2<-testroc(seu=seu2,sce=sce2,outname = "train with createdDB",add=add2)
DBboost<-DBboostTrain(mattest=mattrain2,mfinal=mfinal,ifadd=T)

#没有add做overkill的情况-------
load("seu2_nooverkill_kmeansclustered.robj")
seu2$true<-seu2$label_scds
mattrain<-seu2@meta.data[,c("true","bcds_s","cxds_s","scran_s","dbf_s")]
solo<-read.table("solo.txt")
ddc<-read.table("doubletdetection.txt")
scr<-read.table("scrublet.txt")
add<-cbind(ddc,scr,solo)
colnames(add)<-c("doubletdetection","scrublet","solo")
rownames(add)<-colnames(seu2)
mattrain<-cbind(mattrain,add)
#
DBboostPre(DBboost=DBboost,mattest=mattrain,outname=mfinal)

##整体统计最终各种信息-----
finalmat<-read.csv(paste0("./finalScore.",mfinal,".csv"))
mean(finalmat$chord)
mean(finalmat$chord[(finalmat$true=="Doublet")])
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/aucpr.R")

