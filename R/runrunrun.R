#' Chord
#'
#' remove doublet with scds,bcds,DoubletFinder,doubletCells.
#'
#' @param seu the input seurat object
#' @param seu the input sce object
#' @param seed an integer, random seed
#' @param k an integer,k-means param k
#' @param overkill if True,use overkill
#' @param overkillrate an integer,remove the top ?% doublet-liked cells of any methods' results.(0-1)
#' @param outname The prefix of the output file
#' @param addmethods2 the table merged with other method's scores2
#' @param addmethods1 the table merged with other method's scores1
#' @param mfinal an integer, the number of iterations for which boosting is run or the number of trees to use. Defaults to mfinal=40 iterations.
#' @import Seurat
#' @import scds
#' @import scater
#' @import rsvd
#' @import Rtsne
#' @import cowplot
#' @import DoubletFinder
#' @import scran
#' @import adabag
#' @export
#' @examples chord<-function(seu=NA,doubletrate=NA,mfinal=40,k=20,overkill=T,overkillrate=1,outname="out",seed=1)

#Chord------
chord<-function(
  seu=NA,
  sce=NA,
  doubletrate=NA,
  mfinal=40,
  k=20,
  overkill=T,
  overkillrate=1,
  outname="out",
  seed=1,
  addmethods1=NA,
  addmethods2=NA
  ){

  require(Seurat)
  require(scds)
  require(scater)
  require(rsvd)
  require(Rtsne)
  require(cowplot)
  require(DoubletFinder)
  require(scran)
  require(adabag)

  if (!(is.na(addmethods1)&is.na(addmethods2))) {

    addmethods2<-read.csv(addmethods2,row.names = 1)
    addmethods1<-read.csv(addmethods1,row.names = 1)
    DBboost<-DBboostTrain(mattest=addmethods2,mfinal=mfinal)
    mattestout<-DBboostPre2(DBboost=DBboost,mattest=addmethods1,seu=seu,sce=sce,outname=paste0(outname,mfinal))

    #seu$chord<-mattestout$chord
    #seu$bcds_s<-mattestout$bcds_s
    #seu$cxds_s<-mattestout$cxds_s
    #seu$dbf_s<-mattestout$dbf_s
    #seu$scran_s<-mattestout$scran_s
    #pdf(paste0(outname,"score.pdf"))
    #print(FeaturePlot(seu,features = c("bcds_s","cxds_s","dbf_s","scran_s","chord")))
    #dev.off()

    write.csv(mattestout,file=paste0(outname,"real_score.csv"))
    d<-rownames(mattestout)[order(mattestout$chord,decreasing = T)[1:round(doubletrate*ncol(seu))]]
    write.csv(d,file=paste0(outname,"_doublet.csv"))
    return(d)
  }

  if (is.na(doubletrate)) {
    stop("You need to specify the percentage of cells you want to remove.（0<doubletrate<1）
         for example~ 0.9% per 1000 cells (10X)")
  }
  if (is.na(seu)) {
    stop("You need to input the seurat object")
  }

  set.seed(seed)
  seu@meta.data<-seu@meta.data[,-grep("pANN",colnames(seu@meta.data))]  #Avoid errors resulting from previous results 2021.06.22
  #doubletrate=0.009*ncol(seu)/1000
  sce<-creatSCE(seu=seu)
  sce<-scds(sce=sce)
  seu<-scranDB(seu=seu)
  seu<-DBF(seu=seu,ground_truth = F,doubletrate=doubletrate)

  mattrain<-testroc(seu=seu,sce=sce,outname = "train")
  write.csv(mattrain,file = "real_data.scores.csv")

  seu2<-overkillDB2(seu=seu,sce=sce,doubletrate=doubletrate,seed=seed,k=k,overkill=overkill,overkillrate=overkillrate)
  doubletrate2=sum(seu2$label_scds=="Doublet")/ncol(seu2)
  sce2<-creatSCE(seu=seu2)
  sce2<-scds(sce=sce2)
  seu2<-scranDB(seu=seu2)
  seu2<-DBF(seu=seu2,ground_truth = F,doubletrate=doubletrate2)
  mattrain2<-testroc2(seu=seu2,sce=sce2,outname = "train with createdDB")
  write.csv(mattrain2,file = "simulated_data.scores.csv")

  DBboost<-DBboostTrain(mattest=mattrain2,mfinal=mfinal)
  mattestout<-DBboostPre(DBboost=DBboost,mattest=NA,seu=seu,sce=sce,outname=paste0(outname,mfinal))

  seu$chord<-mattestout$chord
  seu$bcds_s<-mattestout$bcds_s
  seu$cxds_s<-mattestout$cxds_s
  seu$dbf_s<-mattestout$dbf_s
  seu$scran_s<-mattestout$scran_s
  pdf(paste0(outname,"score.pdf"))
  print(FeaturePlot(seu,features = c("bcds_s","cxds_s","dbf_s","scran_s","chord")))
  dev.off()

  write.csv(mattestout,file=paste0(outname,"real_score.csv"))
  d<-rownames(mattestout)[order(mattestout$chord,decreasing = T)[1:round(doubletrate*ncol(seu))]]
  write.csv(d,file=paste0(outname,"_doublet.csv"))
  save(seu,file="seu.robj")
  save(sce,file="sce.robj")
  return(d)
}



