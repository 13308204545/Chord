
#对于真实数据集，待修改
#ggplot aucroc aucpr------------

rocplot<- function(pred, truth, ...){
  
  #打印AUc
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
testroc<-function(seu,sce,outname,add=NULL){
  library(pROC)
  library(ROCR)
  true<-as.vector(sce[,sce$label_scds!="Negative"]$label_scds)
  bcds_s<-sce[,sce$label_scds!="Negative"]$bcds_score
  cxds_s<-sce[,sce$label_scds!="Negative"]$cxds_score
  dbf_s<-seu[,seu$label_scds!="Negative"]@meta.data[,ncol(seu@meta.data)-1]
  scran_s<-seu[,seu$label_scds!="Negative"]@meta.data[,"scran"]
  mattest<-as.data.frame(cbind(true,bcds_s,cxds_s,dbf_s,scran_s))
  mattest$bcds_s<-as.numeric(mattest$bcds_s)
  mattest$cxds_s<-as.numeric(mattest$cxds_s)
  mattest$dbf_s<-as.numeric(mattest$dbf_s)
  mattest$scran_s<-as.numeric(mattest$scran_s)
  #add
  if (!is.null(add)) {
    mattest<-cbind(mattest,add)
  }
  #
  str(mattest)
  
  #roc  aucpr
#  bcds_roc<-roc(formula = true ~ bcds_s, data = mattest)
#  cxds_roc<-roc(formula = true ~ cxds_s, data = mattest)
#  dbf_roc<-roc(formula = true ~ dbf_s, data = mattest)
#  scran_roc<-roc(formula = true ~ scran_s, data = mattest)
  
#  if (!is.null(add)) {
#    for (i in colnames(add)) {
#      roc(formula(paste0("true ~",i)), data = mattest)$auc
#    }
    
#  }
  
  
#  auroc<-c(bcds_roc$auc,cxds_roc$auc,dbf_roc$auc,scran_roc$auc)
  
#  bcds_pr<- prediction(mattest$bcds_s, mattest$true)
#  bcds.aucpr<- as.numeric(performance(bcds_pr, measure = 'aucpr')@y.values)
#  cxds_pr<- prediction(mattest$cxds_s, mattest$true)
#  cxds.aucpr<- as.numeric(performance(cxds_pr, measure = 'aucpr')@y.values)
#  dbf_pr<- prediction(mattest$dbf_s, mattest$true)
#  dbf.aucpr<- as.numeric(performance(dbf_pr, measure = 'aucpr')@y.values)
#  scran_pr<- prediction(mattest$scran_s, mattest$true)
#  scran.aucpr<- as.numeric(performance(scran_pr, measure = 'aucpr')@y.values)
#  aucpr<-c(bcds.aucpr,cxds.aucpr,dbf.aucpr,scran.aucpr)
  
#  write.csv(data.frame(auroc,aucpr,row.names = c("bcds","cxds","dbf","scran")),file = paste0("auc.pr.",outname,".csv"))
  ##在HTOdemux数据集中，bcds和cxds的灵敏度表现较dbf更弱（会漏掉更多TP），但是specificity更准确（FP更低）
  
  #
  #pdf(paste0(outname,"_bcds_cxds_dbf_scran_roc.pdf"))
  #plot.roc(bcds_roc, xlim = c(1,0), ylim = c(0,1), col = "blue", lwd = 2,print.thres = T,print.thres.col = "red",print.auc = T)
  #plot.roc(cxds_roc, xlim = c(1,0), ylim = c(0,1), col = "blue", lwd = 2,print.thres = T,print.thres.col = "red",print.auc = T)
  #plot.roc(dbf_roc, xlim = c(1,0), ylim = c(0,1),col = "blue", lwd = 2,print.thres = T,print.thres.col = "red",print.auc = T)
  #plot.roc(scran_roc, xlim = c(1,0), ylim = c(0,1),col = "blue", lwd = 2,print.thres = T,print.thres.col = "red",print.auc = T)
  #dev.off()
  #cor(mattest[,2:5])

  return(mattest)

}

