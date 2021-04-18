#' overkill the real doublets
#'
#' simulation training set is generated from quality singlet data after filter any doublets detected by any method
#'
#' @param seu the seurat object
#' @param sce the SingleCellExperiment object
#' @param seed random seed
#' @param k k-means param k
#' @param overkill weather use overkill
#' @param overkillrate remove the top ?% doublet-liked cells of any methods' results.

overkillDB2<-function(seu,sce,doubletrate,seed=1,k=20,overkill=T,overkillrate=1){
  require(Seurat)
  require(scds)
  require(scater)
  require(rsvd)
  require(Rtsne)
  require(cowplot)
  seu$bcds_s<-sce$bcds_score
  seu$cxds_s<-sce$cxds_score
  seu$scran_s<-seu@meta.data[,"scran"]
  seu$dbf_s<-seu@meta.data[,grep("pANN",colnames(seu@meta.data))]
  x<-nrow(seu@meta.data)*doubletrate

  seu@meta.data$overkill<-"Singlet"
  seu@meta.data[names(sort(x=seu$bcds_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$cxds_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$scran_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$dbf_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  overkillrate<-sum(seu$overkill%in%"Doublet")/x
  seu2<-seu


  ##PCA+kmeans,delete overkilled doublet
  pcmat<-seu2@reductions$pca@cell.embeddings[,1:30] #first 30pc kmeans
  cl<- kmeans(pcmat,k)   ##k
  seu2$kcluster<-cl$cluster

  pdf("kcluster.pdf")
  print(DimPlot(seu2,group.by = "kcluster"))
  dev.off()


  ##add doublet    #####待优化为针对稀疏矩阵的做,如何最高效率进行随机抽取
  nrow(seu2@meta.data)
  if (overkill==T) {
    mat<-seu2@meta.data[seu2$overkill=="Singlet",]
  }else{
    mat<-seu2@meta.data
  }

  a<-c()
  b<-c()
  x<-c()
  y<-c()
  for (i in 1:k) {
    a<-c(a,sample(x=rownames(mat[mat$kcluster==i,]),size = (doubletrate/(1-doubletrate))*nrow(mat[mat$kcluster==i,]),replace = F))
    b<-c(b,sample(x=rownames(mat[mat$kcluster!=i,]),size = (doubletrate/(1-doubletrate))*nrow(mat[mat$kcluster==i,]),replace = F))
  }  #这里用a=doubletrate/(1-doubletrate)生成双胞数，可以使得原数据双胞率不变
  for (i in 1:length(a)) {
    x<-c(x,rnorm(1,mean=1,sd=0.1))
    y<-c(y,rnorm(1,mean=1,sd=0.1))
  }
  counts<-as.matrix(seu2@assays$RNA@counts)
  counts<-cbind(counts,(counts[,a]*x+counts[,b]*y)/(x+y))##引入高斯分布合

  #----
  seu3<- CreateSeuratObject(counts =counts)
  seu3$label_scds<-"Doublet"   #上标签
  seu3$label_scds[1:nrow(seu2@meta.data)]<-"Singlet" #暂时默认没合成的都是singlet，最后验证会提取true lab标记真实双胞

  seu3@meta.data$overkill<-"Singlet"   ##这里的singlet包含了合成的doublet,方便提取而已
  seu3$overkill[rownames(seu2@meta.data)[seu2$overkill=="Doublet"]]<-"Doublet"

  #-------------

  seu5<-seu3[,seu3$overkill=="Singlet"]
  ##生成异cluster双胞并添加标记后分组------
  overkilled<-seu5
  save(overkilled,file = paste0("overkilled.robj"))
  return(seu5)


}

