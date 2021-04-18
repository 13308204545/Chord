DBF<-function(seu,ground_truth=T,doubletrate){
  library(DoubletFinder)
  library(Seurat)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:10)
  seu<- FindNeighbors(seu, dims = 1:10)
  seu<- FindClusters(seu,verbose = FALSE,resolution=0.4)
  #
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu@meta.data$ClusteringResults
  
  nExp_poi <- round((doubletrate)*length(colnames(seu)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  if(ground_truth==F){
    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(seu, PCs = 1:10, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)
    mpK<-as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  }
  if(ground_truth==T){
    ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
    sweep.res.list_kidney <- paramSweep_v3(seu, PCs = 1:10, sct = FALSE)
    #gt.calls <- seu@meta.data[rownames(sweep.res.list_kidney[[1]]), "label_scds"]  ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
    #sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
    #bcmvn_kidney <- find.pK(sweep.stats_kidney)
  }
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  #seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_150", sct = FALSE)
  return(seu)
}


#seutrain<-DBF(seu=traindata,ground_truth = T)
#seutest<-DBF(seu=testdata,ground_truth = T)





