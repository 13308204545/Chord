##add some packages and test

#install.packages("xgboost")
#library(xgboost)
#install.packages('gbm')
#library(gbm)
#install.packages("lightgbm", repos = "https://cran.r-project.org")
#library(lightgbm)
#mattrain2<-read.csv("simulated_data.scores.csv",row.names = 1)
#mattrain<-read.csv("real_data.scores.csv",row.names = 1)
#mfinal=40
#load("seu.robj")
#load("sce.robj")
#mattest<-mattrain2 #train
#mattest<-mattrain #prediction
#outname="out_gbm"

#####======================run=================================================

##V3:设置了可选择集成方法
DBboostTrain<-function(mattest,mfinal,method="gbm"){
  if (method=="adaboost") {
    require(adabag)
    str(mattest)
    mattest$true<-as.factor(mattest$true)
    DBboost<- boosting(true~.,data = mattest,boos=TRUE, mfinal=mfinal)
    return(DBboost)
  }
  if (method=="xgboost") {
    require(xgboost)
    require(Matrix)
    str(mattest)
    mattest$true<-as.factor(mattest$true)
    train_matrix <- sparse.model.matrix(true ~ .-1, data = mattest)
    train_label <- as.numeric(mattest$true=="Doublet")  # 1:Doublet 0:Singlet
    train_fin <- list(data=train_matrix,label=train_label)
    dtrain<- xgb.DMatrix(data = train_fin$data,label = train_fin$label)
    DBboost <- xgboost(data = dtrain,max_depth=6, eta=0.5, objective='binary:logistic', nround=40)
    importance <- xgb.importance(train_matrix@Dimnames[[2]], model = DBboost)
    head(importance)
    return(DBboost)
  }
  if (method=="lightgbm") {
    require(lightgbm)
    require(Matrix)
    str(mattest)
    mattest$true<-as.factor(mattest$true)
    train_matrix <- sparse.model.matrix(true ~ .-1, data = mattest)
    train_label <- as.numeric(mattest$true=="Doublet")  # 1:Doublet 0:Singlet
    train_fin <- list(data=train_matrix,label=train_label)
    dtrain<- lgb.Dataset(data = train_fin$data,label = train_fin$label)
    train_params <- list(
      num_leaves = 4L
      , learning_rate =0.1
      , objective = "binary"
      , nthread = 2L
      , metric= "l2"
    )
    DBboost <- lightgbm(
      data = train_fin$data
      , params = train_params
      , label = (train_fin$label)
      , nrounds = 2L
    )
    return(DBboost)
  }
  if (method=="gbm") {
    require(gbm)
    str(mattest)
    mattest$true<-as.factor(mattest$true)
    DBboost <- gbm(true~.,distribution = 'gaussian',data=mattest,n.trees=1000,shrinkage = 0.01,cv.folds = 5)
    return(DBboost)
  }
}
DBboostPre<-function(DBboost,seu=seu,sce=sce,mattest=NULL,outname=outname,method="gbm"){

  if (is.null(mattest)) {
    bcds_s<-sce$bcds_score
    cxds_s<-sce$cxds_score
    dbf_s<-seu@meta.data[,grep("pANN",colnames(seu@meta.data))]
#   scran_s<-seu@meta.data[,"scran"]   ######
    mattest<-as.data.frame(cbind(bcds_s,cxds_s,dbf_s))
    mattest$bcds_s<-as.numeric(mattest$bcds_s)
    mattest$cxds_s<-as.numeric(mattest$cxds_s)
    mattest$dbf_s<-as.numeric(mattest$dbf_s)
   #mattest$scran_s<-as.numeric(mattest$scran_s)
    str(mattest)
  }
  if (method=="adaboost"){
    require(adabag)
    pre_DB <- predict(DBboost,newdata = mattest)
    mattest$chord<-pre_DB$prob[,1]-pre_DB$prob[,2]
    write.csv(mattest,file = paste0("finalScore.",outname,".csv"))
    pdf(paste0(outname,"hist.pdf"))
    p<-hist(x = mattest$chord ,breaks = 50, xlim = c(-1,1))
    dev.off()
  }
  if (method=="xgboost"){
    require(xgboost)
    require(Matrix)
    mattest$true=1#凑格式
    test_matrix <- sparse.model.matrix(true ~ .-1, data = mattest)
    test_fin <- list(data=test_matrix,label=mattest$true)
    pre_DB = predict(DBboost,test_fin$data)
    mattest<-mattest[,colnames(mattest)!="true"]
    mattest$chord<-pre_DB
    write.csv(mattest,file = paste0("finalScore.",outname,".csv"))
    pdf(paste0(outname,"hist.pdf"))
    p<-hist(x = mattest$chord ,breaks = 50)
    dev.off()
  }
  if (method=="lightgbm"){
    require(lightgbm)
    require(Matrix)
    mattest$true=1#凑格式
    test_matrix <- sparse.model.matrix(true ~ .-1, data = mattest)
    test_fin <- list(data=test_matrix,label=mattest$true)
    pre_DB = predict(DBboost,test_fin$data)
    mattest<-mattest[,colnames(mattest)!="true"]
    mattest$chord<-pre_DB
    write.csv(mattest,file = paste0("finalScore.",outname,".csv"))
    pdf(paste0(outname,"hist.pdf"))
    p<-hist(x = mattest$chord ,breaks = 50)
    dev.off()
  }
  if (method=="gbm"){
    require(gbm)
    pred<-predict(DBboost,mattest,gbm.perf(DBboost))
    mattest$chord<--pred
    write.csv(mattest,file = paste0("finalScore.",outname,".csv"))
    pdf(paste0(outname,"hist.pdf"))
    p<-hist(x = mattest$chord ,breaks = 50)
    dev.off()
    }
  return(mattest)
}
