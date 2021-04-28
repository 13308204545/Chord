#' translate seurat obj to SCE
#'
#' simulation training set is generated from quality singlet data after filter any doublets detected by any method
#'
#' @param seu the seurat object
creatSCE<-function(seu){
  require(SingleCellExperiment)
  #load(wd)
  # input matrix
  counts_data <- as.matrix(seu@assays$RNA@data)
  col<-seu@meta.data
  sce <- SingleCellExperiment(assays = counts_data,colData = col)
  return(sce)
}


