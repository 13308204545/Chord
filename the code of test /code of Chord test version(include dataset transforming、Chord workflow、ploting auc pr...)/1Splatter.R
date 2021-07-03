library(splatter)
library(scater)
library(Seurat)
set.seed(1)
#sce <- mockSCE()#包自带测试数据
#load('/Users/xiongkexu/Desktop/marked/CNV/remove/reshigh/remove2/Harmony2.robj')
#NPC<-NPC[,sample(length(colnames(NPC)),size = 5000,replace = F)]
#NPC<-NPC[sample(length(rownames(NPC)),size = 500,replace = F),]
## NOTE: Library sizes have been found to be normally distributed instead of log-normal. You may want to check this is correct.
#counts_data <- as.matrix(NPC@assays$RNA@counts)
#col<-NPC@meta.data
#sce<- SingleCellExperiment(assays = list(counts = counts_data),colData = col)
#params <- splatEstimate(sce)
#params

#sim.groups <- splatSimulate(params,group.prob = c(0.3,0.2,0.4,0.1), method = "groups",verbose = FALSE)
#sim.groups <- logNormCounts(sim.groups)
#sim.groups <- runPCA(sim.groups)
#plotPCA(sim.groups, colour_by = "Group")

#test
sim.groups <- splatSimulate(group.prob = c(0.2,0.3,0.4,0.1), method = "groups",verbose = FALSE,nGenes = 8000,batchCells=2000)
sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups)
plotPCA(sim.groups, colour_by = "Group")

setwd("/Users/xiongkexu/Desktop/双胞项目/data/splatter1/")
save(sim.groups,file="splatter1.robj")
counts<-sim.groups@assays@data$TrueCounts

dou<-counts[,sample(1:2000,1)]+counts[,sample(1:2000,1)]/2
for (i in 1:99) {
  dou<-cbind(dou,counts[,sample(1:2000,1)]+counts[,sample(1:2000,1)]/2)
}
counts2<-cbind(counts,dou)
dim(counts2)
seu <- CreateSeuratObject(counts =counts2)
seu$label_scds<-c(rep("Single",2000),rep("Doublet",200))
save(seu,file="splatter2_seu.robj")

setwd("/Users/xiongkexu/Desktop/双胞项目/data/splatter1/")
load("splatter2_seu.robj")
set.seed(1)
seu_test_umi<-sample(colnames(seu),size = round(length(colnames(seu))*(30/100)))
seutest<-seu[,seu_test_umi]
testdata<-seu[,!colnames(seu)%in%seu_test_umi]
set.seed(1)
seu_test_umi<-sample(colnames(seu),size = round(length(colnames(seu))*(70/100)))
seutest<-seu[,seu_test_umi]
traindata<-seu[,colnames(seu)%in%seu_test_umi]



#test2
load(file="/Users/xiongkexu/Desktop/双胞项目/data/splatter1（10%）/splatter1.robj")
setwd("/Users/xiongkexu/Desktop/双胞项目/data/splatter2  (5%)/")

counts<-sim.groups@assays@data$TrueCounts

dou<-counts[,sample(1:2000,1)]+counts[,sample(1:2000,1)]/2
for (i in 1:99) {
  dou<-cbind(dou,counts[,sample(1:2000,1)]+counts[,sample(1:2000,1)]/2)
}
counts2<-cbind(counts,dou)
dim(counts2)
seu <- CreateSeuratObject(counts =counts2)
seu$label_scds<-c(rep("Single",2000),rep("Doublet",100))
save(seu,file="splatter2_seu.robj")

setwd("/Users/xiongkexu/Desktop/双胞项目/data/splatter2  (5%)/")
load("splatter2_seu.robj")
set.seed(1)
seu_test_umi<-sample(colnames(seu),size = round(length(colnames(seu))*(30/100)))
seutest<-seu[,seu_test_umi]
testdata<-seu[,!colnames(seu)%in%seu_test_umi]
set.seed(1)
seu_test_umi<-sample(colnames(seu),size = round(length(colnames(seu))*(70/100)))
seutest<-seu[,seu_test_umi]
traindata<-seu[,colnames(seu)%in%seu_test_umi]

