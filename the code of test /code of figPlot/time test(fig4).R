setwd("A:/双胞项目/软件脚本v3")
source("./1HTOdemux.R")
source("./2SingleCellExperiment_INPUT.R")
source("./3cxds boosting.R")
source("./3DBLFinder.R")
source("./4testauc.R")
source("./5adboosting.R")
source("./1overkillDB.R")
source("./3scran.R")
setwd("G:/工作文件/双胞测试数据/timetest")
load("seu2_nooverkill_kmeansclustered.robj")
library(reticulate)
py_available()
Sys.setenv(RETICULATE_PYTHON="B:/Anaconda3/python.exe")
use_condaenv("B:/Anaconda3/")
use_python("B:/Anaconda3/python.exe")
Sys.which("python")
use_condaenv("base")
py_config()


testall<-function(seu.test){
  timelist<-list()
  
  timelist[["Chord"]]<-system.time(
    {
      doubletrate=0.1
      seu<-seu.test
      sce<-creatSCE(seu=seu)
      sce<-scds(sce=sce)
      seu<-scranDB(seu=seu)
      seu<-DBF(seu=seu,ground_truth = T,doubletrate=doubletrate)
      mattrain<-testroc(seu=seu,sce=sce,outname = "train")
      seu2<-overkillDB2(seu=seu,sce=sce,doubletrate=doubletrate,seed=1,out="all",k=15,overkill=T)
      doubletrate2=sum(seu2$label_scds=="Doublet")/ncol(seu2)
      sce2<-creatSCE(seu=seu2)
      sce2<-scds(sce=sce2)
      seu2<-scranDB(seu=seu2)
      seu2<-DBF(seu=seu2,ground_truth = T,doubletrate=doubletrate2)
      mattrain2<-testroc(seu=seu2,sce=sce2,outname = "train with createdDB")
      DBboost<-DBboostTrain(mattest=mattrain2,mfinal=40)
      DBboostPre(DBboost=DBboost,mattest=mattrain,outname="40")
    }
  )
  library(scds)
  library(scater)
  library(rsvd)
  library(Rtsne)
  library(cowplot)
  timelist[["cxds"]]<-system.time({
    cxds(sce)
  })
  timelist[["bcds"]]<-system.time({
    bcds(sce)
  })
  timelist[["doubletCells"]]<-system.time({
    scranDB(seu=seu)
  })
  timelist[["DoubletFinder"]]<-system.time({
    DBF(seu=seu,ground_truth = T,doubletrate=doubletrate)
  })

  
  doubletdetection<-import('doubletdetection')
  np<-import('numpy')
  timelist[["DoubletDetection"]]<-system.time({
    clf<-doubletdetection$BoostClassifier(n_iters=50L,use_phenograph=F,standard_scaling=T)
    clf$fit(t(as.matrix(seu@assays$RNA@counts)))
  })
  
  import('os')
  scr<-import('scrublet')
  timelist[["scrublet"]]<-system.time({
    scrub=scr$Scrublet(t(as.matrix(seu@assays$RNA@counts)))
    doublet_scores=scrub$scrub_doublets()
    scrub$calculate_doublet_scores()
  })
  
  library(loomR)
  write.csv(as.matrix(seu.test@assays$RNA@counts),file = "count.csv")
  as.loom(x=seu.test,filename = 'seu22.loom',verbose = F,overwrite = T)
  
  timelist[["Solo"]]<-system.time({
    system("solo G:/工作文件/双胞测试数据/HTO8gra/solo_params_example.json G:/工作文件/双胞测试数据/timetest/seu22.loom")
  })
   return(timelist)
}

seu.test<-seu2[,1:1000]
t1000<-testall(seu.test = seu.test)
save(t1000,file = "t1000.robj")
seu.test<-seu2[,1:2000]
t2000<-testall(seu.test = seu.test)
save(t2000,file = "t2000.robj")
seu.test<-seu2[,1:3000]
t3000<-testall(seu.test = seu.test)
save(t3000,file = "t3000.robj")
seu.test<-seu2[,1:4000]
t4000<-testall(seu.test = seu.test)
save(t4000,file = "t4000.robj")
seu.test<-seu2[,1:5000]
t5000<-testall(seu.test = seu.test)
save(t5000,file = "t5000.robj")
seu.test<-seu2[,1:6000]
t6000<-testall(seu.test = seu.test)
save(t5000,file = "t6000.robj")
seu.test<-seu2[,1:7000]
t7000<-testall(seu.test = seu.test)
save(t6000,file = "t7000.robj")
seu.test<-seu2[,1:8000]
t8000<-testall(seu.test = seu.test)
save(t8000,file = "t8000.robj")
seu.test<-seu2[,1:9000]
t9000<-testall(seu.test = seu.test)
save(t9000,file = "t9000.robj")
seu.test<-seu2[,1:10000]
t10000<-testall(seu.test = seu.test)
save(t10000,file = "t10000.robj")
seu.test<-seu2[,1:11000]
t11000<-testall(seu.test = seu.test)
save(t11000,file = "t11000.robj")
seu.test<-seu2[,1:12000]
t12000<-testall(seu.test = seu.test)
save(t12000,file = "t12000.robj")


#整理画图
library(magrittr)
namelist<-Sys.glob("t*.robj")
names(data)
mat<-as.data.frame(matrix(nrow=8,ncol=13))
rownames(mat) =c("Chord","cxds","bcds","doubletCells","DoubletFinder","DoubletDetection","scrublet","Solo" )
colnames(mat)=c((1:12)*1000,"method")

for (i in namelist) {
  x<-load(i)
  data <- eval(parse(text = x))
  mat[,sub("t","",sub(".robj","",i))]<-data.frame(matrix(t(sapply(data,c)),nrow=8))[,3]
}
mat$method<-rownames(mat)
write.csv(mat,file = "timeall.csv")
