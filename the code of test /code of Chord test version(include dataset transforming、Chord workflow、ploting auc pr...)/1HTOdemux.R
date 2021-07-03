#setwd("/Users/xiongkexu/Desktop/双胞项目/data/HTOdemux8/")
#pbmc.umis <- readRDS("./pbmc_umi_mtx.rds")
#pbmc.htos <- readRDS("./pbmc_hto_mtx.rds")

#setwd("/Users/xiongkexu/Desktop/双胞项目/data/HTOdemux12/")
#pbmc.umis <- readRDS("./hto12_umi_mtx.rds")
#pbmc.htos <- readRDS("./hto12_hto_mtx.rds")

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.


HTO8<-function(pbmc.umis,pbmc.htos,seed=1,out="test",type=8){
  library(Seurat)
  if(type==8){
    joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))
    # Subset RNA and HTO counts by joint cell barcodes
    pbmc.umis <- pbmc.umis[, joint.bcs]
    pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])
    # Confirm that the HTO have the correct names
    rownames(pbmc.htos)
    # Setup Seurat object
    pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)
    # Normalize RNA data with log normalization
    pbmc.hashtag <- NormalizeData(pbmc.hashtag)
    # Find and scale variable features
    pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
    pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))
    # Add HTO data as a new assay independent from RNA
    pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
    # Normalize HTO data, here we use centered log-ratio (CLR) transformation
    pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
  }
  if(type==12){
    cells.use <- intersect(rownames(pbmc.htos), colnames(pbmc.umis))
    # Create Seurat object and add HTO data
    pbmc <- CreateSeuratObject(counts = pbmc.umis[, cells.use], min.features = 300)
    pbmc[["HTO"]] <- CreateAssayObject(counts = t(x = pbmc.htos[colnames(pbmc), 1:12]))
    
    # Normalize data
    pbmc <- NormalizeData(pbmc)
    pbmc.hashtag <- NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
    
  }

  
  # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
  # clustering function for large applications You can also play with additional parameters (see
  # documentation for HTODemux()) to adjust the threshold for classification Here we are using the
  # default settings
  pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
  
  # Global classification results
  table(pbmc.hashtag$HTO_classification.global)
  
  # Group cells based on the max HTO signal
  Idents(pbmc.hashtag) <- "HTO_maxID"
  RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:2], ncol = 2)
  
  #outfile
  seu<-pbmc.hashtag
  seu$label_scds<-seu$HTO_classification.global
  seu<-seu[,seu$label_scds!="Negative"]
  
  if(out=="all"){
    return(seu)
  }
  if(out=="train"){
    #splite 30% data to test,70% train
    set.seed(seed)
    seu_test_umi<-sample(colnames(seu),size = round(length(colnames(seu))*(70/100)))
    seutest<-seu[,seu_test_umi]
    seu<-seu[,!colnames(seu)%in%seu_test_umi]
    return(seu)
  }
  if (out=="test") {
    #splite 30% data to test,70% train
    set.seed(seed)
    seu_test_umi<-sample(colnames(seu),size = round(length(colnames(seu))*(30/100)))
    seutest<-seu[,seu_test_umi]
    seu<-seu[,colnames(seu)%in%seu_test_umi]
    return(seu)
  }
}

#testdata<-HTO8(pbmc.umis,pbmc.htos,seed=1,out="test",type=12)
#traindata<-HTO8(pbmc.umis,pbmc.htos,seed=1,out="train",type=12)



