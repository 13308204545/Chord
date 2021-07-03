##V3：保留overkill2函数
# 合成双胞方法改为使用高斯分布rnorm(1,mean=1,sd=0.1)，分别生成两个细胞的比例，再进行加权平均
load("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/colors.robj")
overkillDB2<-function(seu,sce,doubletrate,seed=1,out="all",k=20,overkill=T,overkillrate=1){
  library(Cairo)
  source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/2SingleCellExperiment_INPUT.R")
  #seu sce
  #doubletrate=0.2
  #k=11
  ##delete cell by mat    input:seu sce doubletrat    out:seu2  某些数据集有“NEGATIVE”类型，这里统一删除这些细胞
  seu$bcds_s<-sce[,sce$label_scds!="Negative"]$bcds_score
  seu$cxds_s<-sce[,sce$label_scds!="Negative"]$cxds_score
  seu$scran_s<-seu[,seu$label_scds!="Negative"]@meta.data[,"scran"]
  seu$dbf_s<-seu[,seu$label_scds!="Negative"]@meta.data[,grep("pANN",colnames(seu@meta.data))]
  x<-nrow(seu@meta.data)*doubletrate
  
  seu@meta.data$overkill<-"Singlet"
  seu@meta.data[names(sort(x=seu$bcds_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$cxds_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$scran_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$dbf_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  overkillrate<-sum(seu$overkill%in%"Doublet")/x
  
  #test<-roc(formula = label_scds ~ cxds_s, data = seu@meta.data,plot=T)#看看混淆矩阵
  table(seu@meta.data$overkill,seu$label_scds)
  #require(caret)  #计算混淆矩阵相关信息的包
  #confusionMatrix(as.factor(seu@meta.data[,ncol(seu@meta.data)-4]),as.factor(seu$label_scds),positive = "Doublet")
  #confusionMatrix(as.factor(seu$overkill),as.factor(seu$label_scds),positive = "Doublet")
  
  seu2<-seu
  pdf("score.pdf")
  print(FeaturePlot(seu,features = c("bcds_s","cxds_s","dbf_s","scran_s")))
  print(DimPlot(seu,group.by = "label_scds"))
  dev.off()
  
  pdf("before overkill.pdf")
  print(DimPlot(seu))
  dev.off()
  pdf("after overkill.pdf")
  print(DimPlot(seu2[,seu$overkill=="Singlet"]))
  dev.off()
  
  ##PCA+kmeans,delete overkilled doublet
  pcmat<-seu2@reductions$pca@cell.embeddings[,1:30] #取前30pc做kmeans
  cl<- kmeans(pcmat,k)   ##k值设定很关键
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
  
  seu3$testdoublet<-seu3$label_scds
  seu3$testdoublet[names(seu$label_scds[seu$label_scds=="Doublet"])]<-"doubletold"   #分老双胞和新双胞聚类看分布
  seu3@meta.data$overkill<-"Singlet"   ##这里的singlet包含了合成的doublet,方便提取而已
  seu3$overkill[rownames(seu2@meta.data)[seu2$overkill=="Doublet"]]<-"Doublet"
  save(seu3,file="seu3_DBadded.robj")
  
  #测试合成双胞合理性

  
  seu3<- NormalizeData(seu3, normalization.method = "LogNormalize", scale.factor = 10000)
  seu3<- FindVariableFeatures(seu3, selection.method = "vst", nfeatures = 2000)
  seu3<- ScaleData(seu3)
  seu3<- RunPCA(seu3)
  seu3<- RunUMAP(seu3, dims = 1:20)
  

  
  #seu4<-seu3[,names(seu3$label_scds[seu3$testdoublet%in%c("Doublet","doubletold")])]
  #seu4<-NormalizeData(seu4, normalization.method = "LogNormalize", scale.factor = 10000)
  #seu4<- FindVariableFeatures(seu4, selection.method = "vst", nfeatures = 2000)
  #seu4<- ScaleData(seu4)
  #seu4<- RunPCA(seu4)
  #seu4<- RunUMAP(seu4, dims = 1:20)
  
  CairoPDF("双胞合成效果测试.pdf")
  ####plot1--------
  pt <- 0.6
  alpha <- 0.6
  p <- DimPlot(seu3, reduction='umap', group.by="testdoublet", pt.size=pt) +
    ggtitle("newDoublet") +
    theme_classic(base_size=16) +
    theme(legend.position='n',
          plot.title = element_text(hjust=0.5, size = 22, face='bold'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18, color = 'black'),
          axis.ticks.length=unit(.25, "cm"))
  
  p$layers[[1]]$aes_params$alpha = alpha
  p$layers[[1]]$aes_params$stroke = 0
  print(p)
  
  ####plot2--------
  #pt <- 0.6
  #alpha <- 0.6
  #p <- DimPlot(seu4, reduction='umap', group.by="testdoublet", pt.size=pt) +
  #  ggtitle("newDoublet") +
  #  theme_classic(base_size=16) +
  #  theme(legend.position='n',
  #        plot.title = element_text(hjust=0.5, size = 22, face='bold'),
  #        axis.title = element_text(size = 18),
  #        axis.text = element_text(size = 18, color = 'black'),
  #        axis.ticks.length=unit(.25, "cm"))
  
  #p$layers[[1]]$aes_params$alpha = alpha
  #p$layers[[1]]$aes_params$stroke = 0

  #print(p)
  
  dev.off()
  #-------------

  save(seu2,file = "seu2_nooverkill_kmeansclustered.robj")
  save(seu3,file="seu3_old and new DB.robj")
  
  seu5<-seu3[,seu3$overkill=="Singlet"]
  ##生成异cluster双胞并添加标记后分组------
  set.seed(seed)
  seu_test_umi<-sample(colnames(seu5),size = round(length(colnames(seu5))*(70/100)))
  seutest<-seu5[,seu_test_umi]
  
  if(out=="all"){
    return(seu5)
  }
  
  if(out=="train"){
    #splite 30% data to test,70% train
    seu5<-seu5[,!colnames(seu5)%in%seu_test_umi]
    return(seu5)
  }
  if (out=="test") {
    #splite 30% data to test,70% train
    seu5<-seu5[,colnames(seu5)%in%seu_test_umi]
    return(seu5)
  }
  
}



overkillDB3<-function(seu,sce,doubletrate,seed=1,out="all",k=20,overkill=T,overkillrate=1){
  library(Cairo)
  source("/Users/xiongkexu/Desktop/双胞项目/软件脚本v3/2SingleCellExperiment_INPUT.R")
  #seu sce
  #doubletrate=0.2
  #k=11
  ##delete cell by mat    input:seu sce doubletrat    out:seu2  某些数据集有“NEGATIVE”类型，这里统一删除这些细胞
  seu$bcds_s<-sce[,sce$label_scds!="Negative"]$bcds_score
  seu$cxds_s<-sce[,sce$label_scds!="Negative"]$cxds_score
  seu$scran_s<-seu[,seu$label_scds!="Negative"]@meta.data[,"scran"]
  seu$dbf_s<-seu[,seu$label_scds!="Negative"]@meta.data[,ncol(seu@meta.data)-1]
  x<-nrow(seu@meta.data)*doubletrate
  
  seu@meta.data$overkill<-"Singlet"
  seu@meta.data[names(sort(x=seu$bcds_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$cxds_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$scran_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  seu@meta.data[names(sort(x=seu$dbf_s,decreasing = T))[1:(overkillrate*x)],"overkill"]<-"Doublet"
  overkillrate<-sum(seu$overkill%in%"Doublet")/x
  
  #test<-roc(formula = label_scds ~ cxds_s, data = seu@meta.data,plot=T)#看看混淆矩阵
  table(seu@meta.data$overkill,seu$label_scds)
  require(caret)  #计算混淆矩阵相关信息的包
  #confusionMatrix(as.factor(seu@meta.data[,ncol(seu@meta.data)-4]),as.factor(seu$label_scds),positive = "Doublet")
  confusionMatrix(as.factor(seu$overkill),as.factor(seu$label_scds),positive = "Doublet")
  
  seu2<-seu
  pdf("score.pdf")
  print(FeaturePlot(seu,features = c("bcds_s","cxds_s","dbf_s","scran_s")))
  print(DimPlot(seu,group.by = "label_scds"))
  dev.off()
  
  pdf("before overkill.pdf")
  print(DimPlot(seu))
  dev.off()
  pdf("after overkill.pdf")
  print(DimPlot(seu2[,seu$overkill=="Singlet"]))
  dev.off()
  
  ##PCA+kmeans,delete overkilled doublet
  pcmat<-seu2@reductions$pca@cell.embeddings[,1:30] #取前30pc做kmeans
  cl<- kmeans(pcmat,k)   ##k值设定很关键
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
  matnum<-sample(1:nrow(mat),size = nrow(mat)/2,replace = F)
  matuse<-mat[matnum,]
  matbank<-mat[-matnum,]
  for (i in 1:k) {
    a<-c(a,sample(x=rownames(matuse[matuse$kcluster==i,]),size = (doubletrate/(1-doubletrate))*nrow(matuse[matuse$kcluster==i,]),replace = F))
    b<-c(b,sample(x=rownames(matuse[matuse$kcluster!=i,]),size = (doubletrate/(1-doubletrate))*nrow(matuse[matuse$kcluster==i,]),replace = F))
  }  #这里用a=doubletrate/(1-doubletrate)生成双胞数，可以使得原数据双胞率不变
    for (i in 1:length(a)) {
    x<-c(x,rnorm(1,mean=1,sd=0.1))
    y<-c(y,rnorm(1,mean=1,sd=0.1))
  }
  counts<-as.matrix(seu2@assays$RNA@counts[,!colnames(seu2@assays$RNA@counts)%in%rownames(matuse)])
  countsuse<-as.matrix(seu2@assays$RNA@counts[,!colnames(seu2@assays$RNA@counts)%in%rownames(matbank)])#双胞备选库
  countsuse<-(countsuse[,a]*x+countsuse[,b]*y)/(x+y)##引入高斯分布合
  colnames(countsuse)<-paste0("doub",colnames(countsuse))
  counts<-cbind(counts,countsuse)
  

  
  
  
  #----
  seu3<- CreateSeuratObject(counts =counts)
  seu3$label_scds<-"Doublet"   #上标签
  seu3$label_scds[rownames(seu3@meta.data)%in%rownames(seu2@meta.data)]<-"Singlet" #暂时默认没合成的都是singlet，最后验证会提取true lab标记真实双胞
  ##!!!!3和2这里改过(上面这句)
  seu3$testdoublet<-seu3$label_scds
  seu3$testdoublet[names(seu$label_scds[seu$label_scds=="Doublet"])]<-"Doubletold"   #分老双胞和新双胞聚类看分布
  seu3@meta.data$overkill<-"Singlet"   ##这里的singlet包含了合成的doublet,方便提取而已  
  seu3$overkill[rownames(seu2@meta.data)[seu2$overkill=="Doublet"]]<-"Doublet"
  save(seu3,file="seu3_DBadded.robj")
  
  #测试合成双胞合理性
  
  
  seu3<- NormalizeData(seu3, normalization.method = "LogNormalize", scale.factor = 10000)
  seu3<- FindVariableFeatures(seu3, selection.method = "vst", nfeatures = 2000)
  seu3<- ScaleData(seu3)
  seu3<- RunPCA(seu3)
  seu3<- RunUMAP(seu3, dims = 1:20)
  
  
  
  seu4<-seu3[,names(seu3$label_scds[seu3$testdoublet%in%c("Doublet","Doubletold")])]
  seu4<-NormalizeData(seu4, normalization.method = "LogNormalize", scale.factor = 10000)
  seu4<- FindVariableFeatures(seu4, selection.method = "vst", nfeatures = 2000)
  seu4<- ScaleData(seu4)
  seu4<- RunPCA(seu4)
  seu4<- RunUMAP(seu4, dims = 1:20)
  
  CairoPDF("双胞合成效果测试.pdf")
  ####plot1--------
  pt <- 0.6
  alpha <- 0.6
  p <- DimPlot(seu3, reduction='umap', group.by="testdoublet", pt.size=pt) +
    ggtitle("newDoublet") +
    theme_classic(base_size=16) +
    theme(legend.position='n',
          plot.title = element_text(hjust=0.5, size = 22, face='bold'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18, color = 'black'),
          axis.ticks.length=unit(.25, "cm"))
  
  p$layers[[1]]$aes_params$alpha = alpha
  p$layers[[1]]$aes_params$stroke = 0
  #print(DimPlot(seu3,reduction = "umap",group.by="testdoublet",pt.size = 0.2 ))
  print(p)
  
  ####plot2--------
  pt <- 0.6
  alpha <- 0.6
  p <- DimPlot(seu4, reduction='umap', group.by="testdoublet", pt.size=pt) +
    ggtitle("newDoublet") +
    theme_classic(base_size=16) +
    theme(legend.position='n',
          plot.title = element_text(hjust=0.5, size = 22, face='bold'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18, color = 'black'),
          axis.ticks.length=unit(.25, "cm"))
  
  p$layers[[1]]$aes_params$alpha = alpha
  p$layers[[1]]$aes_params$stroke = 0
  
  print(p)
  
  #print(DimPlot(seu4,reduction = "umap",group.by="testdoublet",pt.size = 0.2 ))
  dev.off()
  #-------------
  
  save(seu2,file = "seu2_nooverkill_kmeansclustered.robj")
  save(seu3,file="seu3_old and new DB.robj")
  
  seu5<-seu3[,seu3$overkill=="Singlet"]
  ##生成异cluster双胞并添加标记后分组------
  set.seed(seed)
  seu_test_umi<-sample(colnames(seu5),size = round(length(colnames(seu5))*(70/100)))
  seutest<-seu5[,seu_test_umi]
  
  if(out=="all"){
    return(seu5)
  }
  
  if(out=="train"){
    #splite 30% data to test,70% train
    seu5<-seu5[,!colnames(seu5)%in%seu_test_umi]
    return(seu5)
  }
  if (out=="test") {
    #splite 30% data to test,70% train
    seu5<-seu5[,colnames(seu5)%in%seu_test_umi]
    return(seu5)
  }
  
}
