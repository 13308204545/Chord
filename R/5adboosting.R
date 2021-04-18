##V2:设置了可选参数
DBboostTrain<-function(mattest,mfinal){
  require(adabag)
  str(mattest)
  mattest$true<-as.factor(mattest$true)
  DBboost<- boosting(true~bcds_s+cxds_s+dbf_s+scran_s,data = mattest,boos=TRUE, mfinal=mfinal)
  return(DBboost)
}
DBboostPre<-function(DBboost,seu=seu,sce=sce,mattest=mattrain,outname="out"){
  require(adabag)

  bcds_s<-sce$bcds_score
  cxds_s<-sce$cxds_score
  dbf_s<-seu@meta.data[,ncol(seu@meta.data)-1]
  scran_s<-seu@meta.data[,"scran"]
  mattest<-as.data.frame(cbind(bcds_s,cxds_s,dbf_s,scran_s))
  mattest$bcds_s<-as.numeric(mattest$bcds_s)
  mattest$cxds_s<-as.numeric(mattest$cxds_s)
  mattest$dbf_s<-as.numeric(mattest$dbf_s)
  mattest$scran_s<-as.numeric(mattest$scran_s)

  str(mattest)
  pre_DB <- predict(DBboost,newdata = mattest)
  mattest$chord<-pre_DB$prob[,1]-pre_DB$prob[,2]

  write.csv(mattest,file = paste0("finalScore.",outname,".csv"))

  pdf(paste0(outname,"hist.pdf"))
  p<-hist(x = mattest$chord ,breaks = 50, xlim = c(-1,1))
  dev.off()
  return(mattest)
}


