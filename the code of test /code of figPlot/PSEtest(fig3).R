#setwd("/Volumes/Transcend/双胞数据/simulated_pse/2/")
setwd("/Volumes/Transcend/双胞数据/simulated_pse/1/")
library(Seurat)
library(monocle)
#手动加载
#seu<-CreateSeuratObject(sim_psudotime_3_sequential[[1]])
#table(sim_psudotime_3_sequential[[2]])
#seu<-seu[,sim_psudotime_3_sequential[[2]]==0]

seu<-CreateSeuratObject(sim_psudotime_bifurcating[[1]])
table(sim_psudotime_bifurcating[[2]])  #这次的数据，生成的双胞也在里面
#seu<-seu[,sim_psudotime_bifurcating[[2]]==0]

seu$label_scds<-"Singlet"
seu$label_scds[which(sim_psudotime_bifurcating[[2]]==1)]<-"Doublet"
#seu$label_scds[which(sim_psudotime_3_sequential[[2]]==1)]<-"Doublet"
#save(seu,file="seuclean.robj")
monocle.run<-function(seu=seu,out=out,color=c("#EE3B3B","lightgray")){
  data <- as(as.matrix(seu@assays$RNA@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = seu@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  monocle_cds <- newCellDataSet(data,phenoData = pd,featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
  HSMM<-monocle_cds
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM,value=F)#
  HSMM <- detectGenes(HSMM, min_expr = 0.1 )
  print(head(HSMM@featureData@data))
  expressed_genes <- row.names(subset(HSMM@featureData@data,num_cells_expressed >=10))
  print(head(pData(HSMM)))
  length(expressed_genes)#查看剩余基
  #HSMM <- setOrderingFilter(HSMM, expressed_genes )
  #plot_ordering_genes(HSMM)
  #采用聚类找差异基因
  HSMM <- reduceDimension(HSMM,
                              max_components = 2,
                              norm_method = 'log',
                              num_dim = 3,
                              reduction_method = 'tSNE',
                              verbose = T)
  HSMM <- clusterCells(HSMM, verbose = F,num_clusters =3)
  #plot_cell_clusters(HSMM, color_by = 'as.factor(Cluster)')
  clustering_DEG_genes <-differentialGeneTest(HSMM[expressed_genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 1)
  HSMM_ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:100]
  HSMM<-setOrderingFilter(HSMM,ordering_genes = HSMM_ordering_genes)
  #plot_ordering_genes(HSMM)
  
  #画图
  HSMM <- reduceDimension(HSMM,method = 'DDRTree')
  HSMM <- orderCells(HSMM)
  pdf(paste0(out,".pdf"))
  print(plot_cell_trajectory(HSMM,color_by = "label_scds",cell_size = 2,cell_link_size=1.5,)+
          scale_color_manual(values = color)+
          theme_bw()+
          theme(legend.position="None",
                panel.grid.major = element_blank(), #主网格线
                panel.grid.minor = element_blank(), #次网格线
                axis.text = element_blank(),
                axis.ticks = element_blank())+
          xlab("")+
          ylab("")
        )
  print(plot_cell_trajectory(HSMM,color_by = "label_scds",cell_size = 2,cell_link_size=1.5,)+
          scale_color_manual(values = color)
  )#提取legend
  dev.off()
  save(HSMM,file=paste0(out,".robj"))
}

#slingone------
slingshot.run<-function(seu=seu,out=out,color=c("#EE3B3B","lightgray")){
  library(SingleCellExperiment)
  library(slingshot)
  library(mclust)
  data <- seu
  count <-seu@assays$RNA@counts
  sce <- SingleCellExperiment(assays = List(counts = count))
  FQnorm <- function(count){ 
    rk <- apply(count,2,rank,ties.method="min")
    count.sort <- apply(count,2,sort)
    refdist <- apply(count.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(count)
    return(norm)
  }
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  embedding.pca <- pca$x[,1:2]
  reducedDims(sce) <- SimpleList(PCA = embedding.pca)
  label.cluster <- Mclust(embedding.pca)$classification
  colData(sce)$GMM <- label.cluster
  sce <- slingshot(sce, clusterLabels = "GMM", reducedDim = "PCA")
  embedding.pca <- as.data.frame(embedding.pca)
  embedding.pca$type <- as.factor(seu$label_scds)
  palette(color)
  pdf(paste0(out,".so.pdf"))
  plot(embedding.pca$PC1, embedding.pca$PC2, col = embedding.pca$type, pch=16, asp = 0,ann = F,xaxt = "n", yaxt ="n",)
  lines(SlingshotDataSet(sce), lwd=4, col="black")
  title(main =out,)
  dev.off()
}

#-----
#monocle.run(seu=seu,out="clean")

#add doublet-------
#
#a<-seu@assays$RNA@counts[,sample(round(0.2*length(seu$orig.ident)))]*rnorm(1,mean=1,sd=0.1)
#b<-seu@assays$RNA@counts[,sample(round(0.2*length(seu$orig.ident)))]*rnorm(1,mean=1,sd=0.1)
#c<-a+b
#d<-cbind(seu@assays$RNA@counts,c)
#rownames(d)<-rownames(seu)
#colnames(d)<-as.character(c(1:ncol(d)))
#seu<-CreateSeuratObject(d)
#seu$label_scds<-"Singlet"
#seu$label_scds[(length(seu$orig.ident)-ncol(c)+1):length(seu$orig.ident)-ncol(c)]<-"Doublet"
#table(seu$label_scds)



#使用原版双胞--------

#V3版run--------
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/1HTOdemux.R")
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/2SingleCellExperiment_INPUT.R")
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/3cxds boosting.R")
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/3DBLFinder.R")
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/4testauc.R")
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/5adboosting.R")
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/1overkillDB.R")
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/3scran.R")

mfinal=40
k=20
overkill=T
overkillrate=1
doubletrate=0.2
system.time(
  { sce<-creatSCE(seu=seu)
  sce<-scds(sce=sce)
  seu<-scranDB(seu=seu)
  seu<-DBF(seu=seu,ground_truth = F,doubletrate=doubletrate)
  
  mattrain<-testroc(seu=seu,sce=sce,outname = "train")#第一次测试,只对有实验标签的使用
  write.csv(mattrain,file="realscore.csv")
  
  #seu2<-overkillDB(seu=seu,sce=sce,doubletrate=doubletrate,seed=1,out="all",k=15)
  seu2<-overkillDB2(seu=seu,sce=sce,doubletrate=doubletrate,seed=1,out="all",k=k,overkill=overkill,overkillrate=overkillrate)
  
  doubletrate2=sum(seu2$label_scds=="Doublet")/ncol(seu2)
  
  sce2<-creatSCE(seu=seu2)
  
  sce2<-scds(sce=sce2)
  
  seu2<-scranDB(seu=seu2)
  
  seu2<-DBF(seu=seu2,ground_truth = F,doubletrate=doubletrate2)
  
  mattrain2<-testroc(seu=seu2,sce=sce2,outname = "train with createdDB")
  
  DBboost<-DBboostTrain(mattest=mattrain2,mfinal=mfinal)
  DBboostPre(DBboost=DBboost,mattest=mattrain,outname=mfinal)
  
  ##整体统计最终各种信息
  finalmat<-read.csv(paste0("./finalScore.",mfinal,".csv"))
  mean(finalmat$chord)
  mean(finalmat$chord[(finalmat$true=="Doublet")])
  source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/aucpr.R")}
)#time:136
save(seu,file="seu.robj")
save(seu2,file="seu2.robj")

library(reticulate)
library(scater)
library(loomR)
library(Seurat)
library(patchwork)
use_condaenv("base")
Sys.setenv(RETICULATE_PYTHON = '/anaconda3/bin/python')
use_python('/anaconda3/bin/python',required = T)

mx<-seu@assays$RNA@counts
write.csv(mx,file = "counts.csv")
#source_python("/Users/xiongkexu/Desktop/python test/3x.py")
system.time(source_python("/Users/xiongkexu/Desktop/python test/1release.py")) #23  这里为了让其有值
system.time(source_python("/Users/xiongkexu/Desktop/python test/2.py")) #2
as.loom(x=seu, filename = "seu22.loom", verbose = FALSE,overwrite = T)

#整理finalmat
finalmat<-read.csv("finalScore.40.csv")
source_python("/Users/xiongkexu/Desktop/python test/4(treat npy)x.py")
solo<-read.table("solo2.txt")
scrublet<-read.table("scrublet.txt")
doubletdetection<-read.table("doubletdetection.txt")
finalmat<-cbind(finalmat[,1:6],doubletdetection,scrublet,solo,finalmat[,7])
colnames(finalmat)<-c("X","true","bcds","cxds","doubletFinder","scran","doubletdetection","scrublet","solo","chord")
#finalmat<-cbind(finalmat[,1:6],doubletdetection,scrublet,finalmat[,7])
#colnames(finalmat)<-c("X","true","bcds","cxds","doubletFinder","scran","doubletdetection","scrublet","chord")
mean(finalmat$chord)
mean(finalmat$chord[(finalmat$true=="Doublet")])
source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/aucpr.R") #这步有把chord变负

#开始剔除
load("seu.robj")
monocle.run(seu=seu,out="contaminated")
slingshot.run(seu=seu,out="contaminated")
seu3<-seu[,order(x=finalmat[,"solo"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]]
monocle.run(seu=seu3,out="solo")
slingshot.run(seu=seu3,out="solo")
seu3<-seu[,order(x=finalmat[,"chord"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]] #先aucpr.R再用这个就改成T
monocle.run(seu=seu3,out="chord")
slingshot.run(seu=seu3,out="chord")
seu3<-subset(seu,subset=label_scds=="Singlet")
monocle.run(seu=seu3,out="clean",color ="lightgray")
slingshot.run(seu=seu3,out="clean",color =c("lightgray","#EE3B3B"))

seu3<-seu[,order(x=finalmat[,"bcds"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]]
monocle.run(seu=seu3,out="bcds")
slingshot.run(seu=seu3,out="bcds")
seu3<-seu[,order(x=finalmat[,"cxds"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]]
monocle.run(seu=seu3,out="cxds")
slingshot.run(seu=seu3,out="cxds")
seu3<-seu[,order(x=finalmat[,"scran"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]]
monocle.run(seu=seu3,out="scran")
slingshot.run(seu=seu3,out="scran")
seu3<-seu[,order(x=finalmat[,"doubletFinder"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]]
monocle.run(seu=seu3,out="doubletFinder")
slingshot.run(seu=seu3,out="doubletFinder")
seu3<-seu[,order(x=finalmat[,"doubletdetection"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]]
monocle.run(seu=seu3,out="doubletdetection")
slingshot.run(seu=seu3,out="doubletdetection")
seu3<-seu[,order(x=finalmat[,"scrublet"],decreasing = T)[round((ncol(seu)*(1/6)+1)):ncol(seu)]]
monocle.run(seu=seu3,out="scrublet")
slingshot.run(seu=seu3,out="scrublet")

table(seu3$label_scds)

#sim_psudotime_bifurcating test doubletdetection

