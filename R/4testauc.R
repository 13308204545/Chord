
#
#ggplot aucroc aucpr------------

rocplot<- function(pred, truth, ...){

  #print AUc
  preauc<-data.frame(pred,truth)
  rocauc<-roc(formula = truth ~ pred, data = preauc)
  predob<- prediction(pred, truth)
  perf.aucpr<- performance(predob, measure = 'aucpr')
  #auroc plot
  perfauc<- performance(predob, 'tpr','fpr')   ## tpr=recall=sensitivity    prec=准确率
  df<- data.frame(x = attributes(perfauc)$x.values[[1]],y = attributes(perfauc)$y.values[[1]])
  p    <- ggplot(data = df)
  p + geom_line(aes(x,y),colour = "yellowgreen",size = 1) +
    geom_ribbon(aes(x,ymin = 0,ymax = y),fill = alpha("yellowgreen",0.5)) +
    labs(title = paste("ROC Curve & AUC:",(round(rocauc$auc,digits = 4)))) +
    xlab("Specificity") +
    ylab("Sensitivity") +
    theme(plot.title = element_text(size = 17))
}

## testroc------------
testroc<-function(seu,sce,outname){
# true<-as.vector(sce[,sce$label_scds!="Negative"]$label_scds)
  bcds_s<-sce$bcds_score
  cxds_s<-sce$cxds_score
  dbf_s<-seu@meta.data[,grep("pANN",colnames(seu@meta.data))]
#  scran_s<-seu@meta.data[,"scran"]   ########
#  mattest<-as.data.frame(cbind(true,bcds_s,cxds_s,dbf_s,scran_s))  #######
  mattest<-as.data.frame(cbind(bcds_s,cxds_s,dbf_s))
  mattest$bcds_s<-as.numeric(mattest$bcds_s)
  mattest$cxds_s<-as.numeric(mattest$cxds_s)
  mattest$dbf_s<-as.numeric(mattest$dbf_s)
#  mattest$scran_s<-as.numeric(mattest$scran_s)  ########
  str(mattest)

  return(mattest)

}

testroc2<-function(seu,sce,outname){
  true<-as.vector(sce$label_scds)
  bcds_s<-sce$bcds_score
  cxds_s<-sce$cxds_score
  dbf_s<-seu@meta.data[,grep("pANN",colnames(seu@meta.data))]
  mattest<-as.data.frame(cbind(true,bcds_s,cxds_s,dbf_s))
  mattest$bcds_s<-as.numeric(mattest$bcds_s)
  mattest$cxds_s<-as.numeric(mattest$cxds_s)
  mattest$dbf_s<-as.numeric(mattest$dbf_s)
  str(mattest)

  return(mattest)

}
