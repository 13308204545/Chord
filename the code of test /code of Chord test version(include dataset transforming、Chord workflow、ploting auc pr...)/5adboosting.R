##V2:设置了可选参数
DBboostTrain<-function(mattest,mfinal,ifadd=F){
  library(adabag)
  str(mattest)
  mattest$true<-as.factor(mattest$true)
  if (ifadd==T) {
    name.method<-""
    for (i in colnames(mattest)[6:ncol(mattest)]) {
      name.method<-paste0(name.method,"+",i)
    }
    DBboost<- boosting(formula(paste0("true ~bcds_s+cxds_s+dbf_s+scran_s",name.method)),data = mattest,boos=TRUE, mfinal=mfinal)
  }else{
    DBboost<- boosting(true~bcds_s+cxds_s+dbf_s+scran_s,data = mattest,boos=TRUE, mfinal=mfinal)
  }
  #预测测试集
  pre_DB <- predict(DBboost,newdata = mattest)
  #将测试集计算所得概率与观测本身取值整合到一起
  obs_DB= data.frame(prob=pre_DB$class,obs=mattest$true)
  #输出混淆矩阵
  table(mattest$true,pre_DB$class,dnn=c("真实值","预测值"))
  #绘制ROC图像
  roc_DBboosting <- roc(mattest$true,(pre_DB$prob[,2]-pre_DB$prob[,1]))
  #pdf("train_roc.pdf")
  #print(plot(roc_DBboosting, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='Adaboost算法ROC曲线'))
  #dev.off()
  return(DBboost)
}
DBboostPre<-function(DBboost,mattest=mattrain,outname="out"){
  library(adabag)
  str(mattest)
  mattest$true<-as.factor(mattest$true)
  #预测测试集
  pre_DB <- predict(DBboost,newdata = mattest)
  #将测试集计算所得概率与观测本身取值整合到一起
  obs_DB= data.frame(prob=pre_DB$class,obs=mattest$true)
  #输出混淆矩阵
  table(mattest$true,pre_DB$class,dnn=c("真实值","预测值"))
  roc_DBboosting <- roc(mattest$true,(pre_DB$prob[,2]-pre_DB$prob[,1])) #利用打分差值作为评分
  mattest$chord<-pre_DB$prob[,2]-pre_DB$prob[,1]
  write.csv(mattest,file = paste0("finalScore.",outname,".csv"))
  pdf(paste0(outname,"vail_roc.pdf"))
  p<-plot(roc_DBboosting, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='Adaboost ROC')
  dev.off()
  return("Finished")
}


#DBboost<-DBboostTrain(mattest=mattrain)
#DBboostPre(DBboost)
