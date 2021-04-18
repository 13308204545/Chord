scranDB<-function(seu=seu){
  require(scran)
  dbl<-doubletCells(seu@assays$RNA@counts)
  seu$scran<-dbl
  return(seu)
}
