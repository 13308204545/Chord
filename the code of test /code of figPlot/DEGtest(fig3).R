#scibet检验 前后差异 
#差异基因 前后差异
#

library(Cairo)
library(ggplot2)
library(Seurat)
library(Matrix)
library(loomR)
library(ggpubr)
library(fpc)
require("RColorBrewer")
display.brewer.all(type = "qual")
color4<-brewer.pal(3, "Set2")
load("/Users/xiongkexu/Desktop/双胞项目/软件脚本v4/color.final.robj")
suppressMessages(library(ROGUE))
suppressMessages(library(tidyverse))
setwd("/Volumes/Transcend/双胞数据/lungCancer2")
#seudata<-Read10X_h5(filename="./NSCLC_EMTAB6149_expression.h5", use.names = TRUE, unique.features = TRUE)
#meta<-read.table("NSCLC_EMTAB6149_CellMetainfo_table.tsv",sep="\t",header = T)
#seu<-CreateSeuratObject(counts=seudata,meta.data = meta)
seu <- connect(filename = "Thienpont_Tumors_52k_v4_R_fixed.loom", mode = "r")
seu <- as.Seurat(seu)
# 拆分保存
#for (i in 1:5) {
#  seusub<-subset(seu,subset = PatientNumber==i)
#  save(seusub,file = paste0("seu",i,".robj"))
#}

#
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunTSNE(seu, dims = 1:30)

markers<-FindAllMarkers(seu,group.by="Clusterings.1")
save(seu,file = "seu.robj",version=2)
load("seu.robj")

x<-read.csv("./1/chord.csv",row.names = 1)
seu$doublet<-0
seu$doublet[rownames(x)[which(x==1)]]<-1
x<-read.csv("./2/chord.csv",row.names = 1)
seu$doublet[rownames(x)[which(x==1)]]<-1
x<-read.csv("./3/chord.csv",row.names = 1)
seu$doublet[rownames(x)[which(x==1)]]<-1
x<-read.csv("./4/chord.csv",row.names = 1)
seu$doublet[rownames(x)[which(x==1)]]<-1
x<-read.csv("./5/chord.csv",row.names = 1)
seu$doublet[rownames(x)[which(x==1)]]<-1

seu$doublet[seu$doublet==1]="doublet"
seu$doublet[seu$doublet==0]="singlet"

seu$seurat_clusters<-seu$Clusterings.1
seu@reductions$tsne@cell.embeddings<-cbind(tSNE_1=seu$Embeddings_X.0,tSNE_2=seu$Embeddings_Y.0)
DimPlot(seu,reduction ="tsne",group.by = "seurat_clusters")

seu$Celltype<-"x"
seu$Celltype[seu$seurat_clusters==7]="Epithelial"
seu$Celltype[seu$seurat_clusters==6]="Alveolar"
seu$Celltype[seu$seurat_clusters==5]="Fibroblasts"
seu$Celltype[seu$seurat_clusters==4]="Cancer"
seu$Celltype[seu$seurat_clusters==3]="Endothelial"
seu$Celltype[seu$seurat_clusters==2]="Myeloid"
seu$Celltype[seu$seurat_clusters==1]="B cells"
seu$Celltype[seu$seurat_clusters==0]="T cells"


color<-as.character(colorlist[[3]])
color2<-c("#EE3B3B","lightgrey")
pdf("Fig4A.pdf",width = 14,height = 14)
alpha=1
p<-DimPlot(seu,group.by = "Celltype",raster = F ,pt.size = 0.6,cols =brewer.pal(n = 8,name = "Set2")[c(7,2,5,4,8,6,3,1)] )+
  ggtitle("Celltype") +
  theme_classic(base_size=16) +
  theme(legend.position='n',
        plot.title = element_text(hjust=0.5, size = 22, face='bold'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = 'black'),
        axis.ticks.length=unit(.25, "cm"))

p$layers[[1]]$aes_params$alpha = alpha
p$layers[[1]]$aes_params$stroke = 0
print(p)

p<-DimPlot(seu,group.by  = "doublet", reduction = "umap",label = F ,pt.size = 0.6,raster = F,cols =color2 ) +
  ggtitle("Doublet") +
  theme_classic(base_size=16) +
  theme(legend.position='n',
        plot.title = element_text(hjust=0.5, size = 22, face='bold'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = 'black'),
        axis.ticks.length=unit(.25, "cm"))

p$layers[[1]]$aes_params$alpha = alpha
p$layers[[1]]$aes_params$stroke = 0
print(p)
dev.off()

DimPlot(seu,group.by = "Celltype",raster = F ,pt.size = 0.6,cols =brewer.pal(n = 8,name = "Set2")[c(7,2,5,4,8,6,3,1)] )


CairoPDF("markerPlot.pdf")
FeaturePlot(seu,features = "CLDN18")
FeaturePlot(seu,features = "CLDN5")
FeaturePlot(seu,features = "CAPS")
FeaturePlot(seu,features = "COL1A1")
FeaturePlot(seu,features = "CD79A")
FeaturePlot(seu,features = "LYZ")
FeaturePlot(seu,features = "CD3D")
FeaturePlot(seu,features = "EPCAM")
dev.off()

#处理,生成删除后的细胞数据
seu.killed<-subset(seu,subset = doublet=="singlet")
seu.table<-table(seu$Celltype)
seu.killed.table<-table(seu.killed$Celltype)
seu.doublet<-subset(seu,subset= doublet=="doublet")
table(seu.doublet$Celltype)/table(seu$Celltype)
filter.table<-rbind(lung=round(seu.table/sum(seu.table),digits = 4),lung.filtered=round(seu.killed.table/sum(seu.killed.table),digits = 4),"doublet/cluster"=round(table(seu.doublet$Celltype)/table(seu$Celltype),digits = 4))
write.csv(filter.table,file="cells.csv")
save(seu.killed,file = "seu.killed.robj",version = 2)
load(file = "seu.killed.robj")

df<-as.data.frame(t(rbind(Frequency=filter.table[3,],Celltype=colnames(filter.table),Count=table(seu.doublet$Celltype))))
df$Frequency<-as.numeric(df$Frequency)
df$Count<-as.numeric(df$Count)
df$Frequency<-paste0(df$Frequency*100,"%")

pdf("双胞比例.pdf",height = 7,width = 7)
df$Celltype<- factor(df$Celltype,levels=c("T cells","B cells","Myeloid","Endothelial","Cancer","Fibroblasts","Alveolar","Epithelial"))
p<-ggplot(df, aes(x=Celltype, y=Count)) + geom_bar(aes(fill=Celltype),position="dodge", stat="identity")+scale_fill_manual(values =brewer.pal(n = 8,name = "Set2"))
p+geom_text(aes(label = Frequency), size = 3, hjust = 0.5, vjust = -1, position ="stack")+theme_classic()
dev.off()

#rogue
rogue.res <- rogue(seu@assays$RNA@counts, labels = seu@meta.data$Celltype, samples =seu@meta.data$PatientNumber, platform = "UMI", span = 0.6)
write.csv(rogue.res,file="rogue(ori).csv")
rogue.res2 <- rogue(seu.killed@assays$RNA@counts, labels = seu.killed@meta.data$Celltype, samples =seu.killed@meta.data$PatientNumber, platform = "UMI", span = 0.6)
write.csv(rogue.res2,file="rogue(filtered).csv")

p<-rogue.boxplot(rogue.res)
p$data
p2<-rogue.boxplot(rogue.res2)
p2$data

p$data$group<-"orignal"
p2$data$group<-"filtered"
pdata<-rbind(p$data,p2$data)

pdata$group<-factor(pdata$group,levels = c("orignal","filtered")) 

pdf("ROGUE.pdf")
pp <- ggboxplot(pdata, x = "clusters", y = "ROGUE",
               color = "group", palette = "jco", 
               add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
pp #+ stat_compare_means(aes(group = group), label = "p.signif",paired=T,method ="wilcox.test")


pp <- ggboxplot(pdata, x = "group", y = "ROGUE",
                color = "group", palette = "jco", 
                add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
pp + stat_compare_means(aes(group = group), label = "p.signif",method ="wilcox.test" ,paired=T)

dev.off()

pdf("rogue.sample.pdf")
for (i in 1:5) {
  p<-rogue.boxplot(rogue.res[i,])
  p$data
  p2<-rogue.boxplot(rogue.res2[i,])
  p2$data
  p$data$group<-"orignal"
  p2$data$group<-"filtered"
  pdata<-rbind(p$data,p2$data)
  pdata$group<-factor(pdata$group,levels = c("orignal","filtered")) 
  
  pp <- ggboxplot(pdata, x = "group", y = "ROGUE",
                  color = "group", palette = "jco", 
                  add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  print(pp + stat_compare_means(aes(group = group), label = "p.signif",method ="wilcox.test" ,paired=T))
}
dev.off()

##用umap出图-------
load(file = "seu.robj")
load(file = "seu.killed.robj")
seu <- RunUMAP(seu, dims = 1:30)
seu.killed <- RunUMAP(seu.killed, dims = 1:30)

pdf("Fig4A_umap.pdf",width = 14,height = 14)
alpha=1
p<-DimPlot(seu,group.by = "Celltype",reduction = "umap",raster = F ,pt.size = 0.6)+
  ggtitle("Celltype") +
  theme_classic(base_size=16) +
  theme(legend.position='n',
        plot.title = element_text(hjust=0.5, size = 22, face='bold'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = 'black'),
        axis.ticks.length=unit(.25, "cm"))

p$layers[[1]]$aes_params$alpha = alpha
p$layers[[1]]$aes_params$stroke = 0
print(p)

p<-DimPlot(seu,group.by  = "doublet", reduction = "umap",label = F ,pt.size = 0.6,raster = F) +
  ggtitle("Doublet") +
  theme_classic(base_size=16) +
  theme(legend.position='n',
        plot.title = element_text(hjust=0.5, size = 22, face='bold'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = 'black'),
        axis.ticks.length=unit(.25, "cm"))

p$layers[[1]]$aes_params$alpha = alpha
p$layers[[1]]$aes_params$stroke = 0
print(p)
dev.off()

save(seu,file = "seu.robj")
save(seu.killed,file = "seu.killed.robj")

#待定，计算单个亚群的rogue
#rogue.res3 <- rogue(seu.killed@assays$RNA@counts, labels = seu.killed@meta.data$Celltype, samples =seu.killed@meta.data$Celltype, platform = "UMI", span = 0.6)
seu@active.ident<-as.factor(seu$Celltype)
markers.seu<-FindAllMarkers(seu,logfc.threshold = 0.25,test.use = "wilcox")
seu.killed@active.ident<-as.factor(seu.killed$Celltype)
markers.seu.killed<-FindAllMarkers(seu.killed,logfc.threshold = 0.25,test.use = "wilcox")
save(markers.seu,file = "marker.seu.robj")
save(markers.seu.killed,file = "marker.seu.killed.robj")

load("marker.seu.robj")
load("marker.seu.killed.robj")
#markers.seu2<-markers.seu[which(markers.seu$p_val_adj<0.05),]
#markers.seu.killed2<-markers.seu.killed[which(markers.seu.killed$p_val_adj<0.05),]
markers.seu2<-markers.seu[which(abs(markers.seu$avg_logFC)>0.693),]
markers.seu.killed2<-markers.seu.killed[which(abs(markers.seu.killed$avg_logFC)>0.693),]
#markers.seu2<-markers.seu
#markers.seu.killed2<-markers.seu.killed


markers.seu2.T<-markers.seu2[which(markers.seu2$cluster=="T cells"),]
markers.seu.killed2.T<-markers.seu.killed2[which(markers.seu.killed2$cluster=="T cells"),]
markers.seu2.B<-markers.seu2[which(markers.seu2$cluster=="B cells"),]
markers.seu.killed2.B<-markers.seu.killed2[which(markers.seu.killed2$cluster=="B cells"),]
markers.seu2.Alveolar<-markers.seu2[which(markers.seu2$cluster=="Alveolar"),]
markers.seu.killed2.Alveolar<-markers.seu.killed2[which(markers.seu.killed2$cluster=="Alveolar"),]
markers.seu2.Alveolar<-markers.seu2[which(markers.seu2$cluster=="Alveolar"),]
markers.seu.killed2.Alveolar<-markers.seu.killed2[which(markers.seu.killed2$cluster=="Alveolar"),]
markers.seu2.Cancer<-markers.seu2[which(markers.seu2$cluster=="Cancer"),]
markers.seu.killed2.Cancer<-markers.seu.killed2[which(markers.seu.killed2$cluster=="Cancer"),]
markers.seu2.Endothelial<-markers.seu2[which(markers.seu2$cluster=="Endothelial"),]
markers.seu.killed2.Endothelial<-markers.seu.killed2[which(markers.seu.killed2$cluster=="Endothelial"),]
markers.seu2.Epithelial<-markers.seu2[which(markers.seu2$cluster=="Epithelial"),]
markers.seu.killed2.Epithelial<-markers.seu.killed2[which(markers.seu.killed2$cluster=="Epithelial"),]
markers.seu2.Fibroblasts<-markers.seu2[which(markers.seu2$cluster=="Fibroblasts"),]
markers.seu.killed2.Fibroblasts<-markers.seu.killed2[which(markers.seu.killed2$cluster=="Fibroblasts"),]
markers.seu2.Myeloid<-markers.seu2[which(markers.seu2$cluster=="Myeloid"),]
markers.seu.killed2.Myeloid<-markers.seu.killed2[which(markers.seu.killed2$cluster=="Myeloid"),]


markers.k.new.Alveolar<-markers.seu.killed2.Alveolar[setdiff(markers.seu.killed2.Alveolar$gene,markers.seu2.Alveolar$gene),]
markers.o.new.Alveolar<-markers.seu2.Alveolar[setdiff(markers.seu2.Alveolar$gene,markers.seu.killed2.Alveolar$gene),]
write.csv(markers.k.new.Alveolar,file = "add.Alveolar.csv")
#k:EML4(NK-T marker) 
markers.k.new.T<-markers.seu.killed2.T[setdiff(markers.seu.killed2.T$gene,markers.seu2.T$gene),]
markers.o.new.T<-markers.seu2.T[setdiff(markers.seu2.T$gene,markers.seu.killed2.T$gene),]
write.csv(markers.k.new.T,file = "add.T.csv")
markers.k.new.B<-markers.seu.killed2.B[setdiff(markers.seu.killed2.B$gene,markers.seu2.B$gene),]
markers.o.new.B<-markers.seu2.B[setdiff(markers.seu2.B$gene,markers.seu.killed2.B$gene),]
write.csv(markers.k.new.B,file = "add.B.csv")
markers.k.new.Cancer<-markers.seu.killed2.Cancer[setdiff(markers.seu.killed2.Cancer$gene,markers.seu2.Cancer$gene),]
markers.o.new.Cancer<-markers.seu2.Cancer[setdiff(markers.seu2.Cancer$gene,markers.seu.killed2.Cancer$gene),]
write.csv(markers.k.new.Cancer,file = "add.Cancer.csv")
markers.k.new.Endothelial<-markers.seu.killed2.Endothelial[setdiff(markers.seu.killed2.Endothelial$gene,markers.seu2.Endothelial$gene),]
markers.o.new.Endothelial<-markers.seu2.Endothelial[setdiff(markers.seu2.Endothelial$gene,markers.seu.killed2.Endothelial$gene),]
write.csv(markers.k.new.Endothelial,file = "add.Endothelial.csv")
markers.k.new.Fibroblasts<-markers.seu.killed2.Fibroblasts[setdiff(markers.seu.killed2.Fibroblasts$gene,markers.seu2.Fibroblasts$gene),]
markers.o.new.Fibroblasts<-markers.seu2.Fibroblasts[setdiff(markers.seu2.Fibroblasts$gene,markers.seu.killed2.Fibroblasts$gene),]
write.csv(markers.k.new.Fibroblasts,file = "add.Fibroblasts.csv")
markers.k.new.Epithelial<-markers.seu.killed2.Epithelial[setdiff(markers.seu.killed2.Epithelial$gene,markers.seu2.Epithelial$gene),]
markers.o.new.Epithelial<-markers.seu2.Epithelial[setdiff(markers.seu2.Epithelial$gene,markers.seu.killed2.Epithelial$gene),]
write.csv(markers.k.new.Epithelial,file = "add.Epithelial.csv")
markers.k.new.Myeloid<-markers.seu.killed2.Myeloid[setdiff(markers.seu.killed2.Myeloid$gene,markers.seu2.Myeloid$gene),]
markers.o.new.Myeloid<-markers.seu2.Myeloid[setdiff(markers.seu2.Myeloid$gene,markers.seu.killed2.Myeloid$gene),]
write.csv(markers.k.new.Myeloid,file = "add.Myeloid.csv")

##正常梯度
df<-data.frame(cutoff=0,ngene=0,stage=0)
for (cutoff in c(0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75)) {
  df<-rbind(df,data.frame(cutoff=cutoff,ngene=nrow(markers.seu[which(abs(markers.seu$avg_logFC)>cutoff),]),stage="nofiltered"))
  df<-rbind(df,data.frame(cutoff=cutoff,ngene=nrow(markers.seu.killed[which(abs(markers.seu.killed$avg_logFC)>cutoff),]),stage="filtered"))
  
}
df<-df[-1,]
p<-ggplot(data = df,aes(x=cutoff,y=ngene,group = stage,color=stage))+   #linetype=linetypek单独设置虚线
  geom_point()+
  geom_line()+
  xlab("cutoff logfc")+#横坐标名称
  ylab("ngene")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
        #legend.position = c(.154,.167),#更改图例的位置，放至图内部的左上角
        legend.box.background = element_rect(color="black"))+#为图例田间边框线
  expand_limits(y=c(0.2,1))+#更改横坐标刻度值(离散型)
  scale_color_manual(values =color4[1:2])+
  scale_x_continuous(breaks=seq(min(df$cutoff),max(df$cutoff),0.1))
pdf("DEGs.pdf")
p
dev.off()

##只取正向的DEGs
df<-data.frame(cutoff=0,ngene=0,stage=0)
for (cutoff in c(0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75)) {
  df<-rbind(df,data.frame(cutoff=cutoff,ngene=nrow(markers.seu[which(markers.seu$avg_logFC>cutoff),]),stage="nofiltered"))
  df<-rbind(df,data.frame(cutoff=cutoff,ngene=nrow(markers.seu.killed[which(markers.seu.killed$avg_logFC>cutoff),]),stage="filtered"))
}
df<-df[-1,]
p<-ggplot(data = df,aes(x=cutoff,y=ngene,group = stage,color=stage))+   #linetype=linetypek单独设置虚线
  geom_point()+
  geom_line()+
  xlab("cutoff logfc")+#横坐标名称
  ylab("ngene")+#纵坐标名称
  theme_bw() +#去掉背景灰色
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
        #legend.position = c(.154,.167),#更改图例的位置，放至图内部的左上角
        legend.box.background = element_rect(color="black"))+#为图例田间边框线
  expand_limits(y=c(0.2,1))+#更改横坐标刻度值(离散型)
  scale_color_manual(values =color4[1:2])+
  scale_x_continuous(breaks=seq(min(df$cutoff),max(df$cutoff),0.1))
pdf("DEGs(positive).pdf")
p
dev.off()

pdf("DEGs.hist.pdf")
hist(markers.seu$avg_logFC,breaks = 200)
hist(markers.seu.killed$avg_logFC,breaks = 200)
dev.off()

#时序
slingshot.run<-function(seu=seu,out=out,color=c("#EE3B3B","lightgray")){
  library(SingleCellExperiment)
  library(slingshot)
  library(mclust)
  data <- seu
  count <-seu@assays$RNA@counts
  sce <- SingleCellExperiment(assays = List(counts = count))
  FQnorm <- function(count){ 
    rk <- apply(count,2,rank,ties.method="min")
    count.sort <- apply(count,2,sort)
    refdist <- apply(count.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(count)
    return(norm)
  }
  assays(sce)$norm <- FQnorm(assays(sce)$counts)
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  embedding.pca <- pca$x[,1:2]
  reducedDims(sce) <- SimpleList(PCA = embedding.pca)
  label.cluster <- Mclust(embedding.pca)$classification
  colData(sce)$GMM <- label.cluster
  sce <- slingshot(sce, clusterLabels = "GMM", reducedDim = "PCA")
  embedding.pca <- as.data.frame(embedding.pca)
  embedding.pca$type <- as.factor(seu$doublet)
  save(sce,file = paste0(out,".robj"))
  palette(color)
  pdf(paste0(out,".so.pdf"))
  plot(embedding.pca$PC1, embedding.pca$PC2, col = embedding.pca$type, pch=16, asp = 0,ann = F,xaxt = "n", yaxt ="n",)
  lines(SlingshotDataSet(sce), lwd=4, col="black")
  title(main =out,)
  dev.off()
}
slingshot.run(seu = seu,out = "lung2",color=c("#EE3B3B","lightgray"))

#亚群处理前后
umidf<-seu@meta.data[,c("nCount_RNA","Celltype","doublet")]
umidf$adj<-umidf$doublet
umidf$adj[umidf$doublet=="doublet"]<-"doublet"
pdf("fig4g.pdf")
pp <- ggboxplot(umidf, x = "Celltype", y = "nCount_RNA",
                color = "doublet",group="doublet", palette = "jco", 
                add = "jitter",add.params = list(size=0.03),outlier.shape = NA)+ # palette可以按照期刊选择相应的配色，如"npg"等
  scale_color_manual(values = c("#EE3B3B","lightgray"))
pp + stat_compare_means(aes(group= doublet), label = "p.signif")
dev.off()
