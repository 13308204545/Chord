#' run bcds cxds
#'
#' run bcds cxds and return SCE object
#'
#' @param sce the input SCE object
scds<-function(sce,
               cxds.ntop=NA,
               cxds.binThresh=NA,
               bcds.ntop=NA,
               bcds.srat=NA){
  require(scds)
  require(scater)
  require(rsvd)
  require(Rtsne)
  require(cowplot)
  sce = cxds(sce,retRes = TRUE,ntop = if (is.na(cxds.ntop)){500}else{cxds.ntop},binThresh =if (is.na(cxds.binThresh)){0}else{cxds.binThresh})
  sce = bcds(sce,retRes = TRUE,verb=TRUE,ntop = if (is.na(bcds.ntop)){500}else{bcds.ntop},srat =  if (is.na(bcds.srat)){1}else{bcds.srat})
  return(sce)
}

##函数--------
cxds_getTopPairs <- function(sce,n=100){
  #=======================================

  ind = rowData(sce)$cxds_hvg_ordr[!is.na(rowData(sce)$cxds_hvg_ordr)]
  Bp  = counts(sce)[ind,] > metadata(sce)$cxds$binThresh
  imp = Matrix::t(Matrix::t(Bp) * sce$cxds_score)
  imp = imp  %*% Matrix::t(Bp)
  imp = imp * metadata(sce)$cxds$S
  imp = as.matrix(imp)
  res = list(imp=imp)

  colnames(imp) = seq_len(ncol(imp))
  rownames(imp) = colnames(imp)
  topPrs = matrix(NA,nrow=n,ncol=2)
  for(i in seq_len(n)){
    pr = which(imp == max(imp) ,arr.ind=TRUE)[c(1,2)]
    topPrs[i,] = as.integer(colnames(imp)[pr])
    imp = imp[-pr,]
    imp = imp[,-pr]
  }
  res$topPairs = topPrs
}

get_dblCalls_ROC <- function(scrs_real, scrs_sim, rel_loss=1){
  #=============================================================

  p    = length(scrs_sim)/(length(scrs_sim)+length(scrs_real))
  rc   = roc(response=c(rep(0,length(scrs_real)),rep(1,length(scrs_sim))),predictor=c(scrs_real,scrs_sim))
  sens = rc$sensitivities
  spec = rc$specificities
  #- use youden to find ROC cutoff
  r       = p/rel_loss/(1-p)
  cut_ind = which.max(sens+r*spec-1)
  thresh  = rc$thresholds[cut_ind]

  ndbl    = sum(scrs_real >= thresh)
  fram    = mean(scrs_sim < thresh) #- fraction of sim doublets missed

  res = c(ndbl,thresh,fram)
  names(res) = c("number","threshold","fnr") #- false negative rate
  return(res)
}
get_dblCalls_dist <- function(scrs_real,scrs_sim, type="balanced"){
  #==================================================================

  #- do "balanced errpr"

  if(type == "balanced"){
    #======================

    es  = ecdf(scrs_sim)
    er  = ecdf(scrs_real)
    rtf = function(thresh) 1-er(thresh) #- right tail; decreases mono with arg
    ltf = function(thresh) es(thresh)   #- left tail; increases mono with arg

    res        = uniroot(function(x)ltf(x) - rtf(x), c(min(scrs_real),max(scrs_real)))
    res_val    = res$root
  } else {
    if(!is.numeric(type)) stop("invalid type argument\n")
    res_val = quantile(scrs_sim,prob=type)
  }
  res_ndl    = sum(scrs_real >= res_val)
  res_fnr    = mean(scrs_sim < res_val)
  res        = c(res_ndl,res_val,res_fnr)
  names(res) = c("number","threshold","fnr") #- false negative rate
  return(res)
}
get_dblCalls_ALL <- function(scrs_real,scrs_sim,rel_loss=1){
  #=========================================================
  est_dbl = rbind( get_dblCalls_ROC( scrs_real,scrs_sim,rel_loss),
                   get_dblCalls_dist(scrs_real,scrs_sim,"balanced"),
                   get_dblCalls_dist(scrs_real,scrs_sim,0.1),
                   get_dblCalls_dist(scrs_real,scrs_sim,0.01))
  rownames(est_dbl)  = c("youden","balanced","0.1","0.01")
  return(est_dbl)
}





