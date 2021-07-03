scranDB<-function(seu=seu){
  library(scran)
  dbl<-doubletCells(seu@assays$RNA@counts)
  seu$scran<-dbl
  return(seu)
}
