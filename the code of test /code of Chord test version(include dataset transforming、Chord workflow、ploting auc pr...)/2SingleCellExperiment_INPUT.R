##SingleCellExperiment输入准备
##输入文件：表达矩阵，一列
#wd="/Users/xiongkexu/Desktop/marked/CNV/remove/reshigh/remove2/Harmony.robj"
#test,input from seurat
creatSCE<-function(seu){
  library(SingleCellExperiment)
  #load(wd)
  # 输入数据的表达矩阵
  counts_data <- seu@assays$RNA@data
  col<-seu@meta.data
  # 构建参考数据SingleCellExperiment对象
  #sce <- SingleCellExperiment(assays = list(counts = counts_data))
  sce <- SingleCellExperiment(assays = list(counts = counts_data),colData = col)  #版本不一样可能需要改为assasyData
  return(sce)
}



##2选1
#scetest<-creatSCE(seu=testdata)
#scetrain<-creatSCE(seu=traindata)
