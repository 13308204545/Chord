# Chord
remove single cells Doublets by integrating tools! 
Chord uses the gbm/AdBoost algorithm to integrate different methods for stable and accurate doublets filtered results. 

parameter:

**method** the boost method ("adaboost" or "gbm")

**seu** the input seurat object

**sce** the input sce object

**seed** an integer, random seed

**k** an integer,k-means param k

**overkill**  if True,use overkill

**overkillrate**  an integer,remove the top ?% doublet-liked cells of any methods' results.(0-1)

**outname** The prefix of the output file

**addmethods2** the table merged with other method's scores2

**addmethods1** the table merged with other method's scores1

**mfinal**  an integer, the number of iterations for which boosting is run or the number of trees to use. Defaults to mfinal=40 iterations.(only works when method="adaboost")

**overkilllist**  a vector of cells to be remove in overkill

**adddoublt** doubletrate of cells to be simulate

**cxds.ntop** integer, indimessageing number of top variance genes to consider. Default: 500

**cxds.binThresh**  integer, minimum counts to consider a gene "present" in a cell. Default: 0

**bcds.ntop** integer, indicating number of top variance genes to consider. Default: 500

**bcds.srat** numeric, indicating ratio between orginal number of "cells" and simulated doublets; Default: 1

**dbf.PCs** Number of statistically-significant principal components (e.g., as estimated from PC elbow plot); Default: 1:10

**dbf.pN**  The number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 0.25, based on observation that DoubletFinder performance is largely pN-invariant (see McGinnis, Murrow and Gartner 2019, Cell Systems).

**dbf.pK**  The PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. Optimal pK values can be determined using mean-variance-normalized bimodality coefficient.


## Install:
```R
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') 

devtools::install_github('kostkalab/scds',ref="master")

install.packages("adabag")

install.packages("gbm")

devtools::install_github("13308204545/Chord") 
   
```
## Quick start:
```R
chord（seu="input seurat object",doubletrate="estimated doubletrate",overkill=T,outname="the name you want"）
```
Q:how to estimate doubletrate? 

A:It depends on the number of cells in the sample. 10X can be referred：doubletrate = ~0.9% per 1,000 cells.  

Q:how to remove doublets 

A:The doublets' barcodes are in the file "outname_doublets.csv" 

## Choose overkill combination
```R
You can choose any combination of methods for overkill by setting the overkilllist parameter (a vector of cells to be remove in overkill)

```
## Boost more methods:
1.Using any method to evaluate the dataset "overkilled.robj", adding the results of socres to "simulated_data.scores.csv".

![image](https://github.com/13308204545/Chord/blob/main/pictures/readme1.png)

2.Using any method to evaluate the dataset "seu.robj", adding the results of socres to "simulated_data.scores.csv".

![image](https://github.com/13308204545/Chord/blob/main/pictures/readme2.png)

3.In the same dir, run the codes:
```R
load("seu.robj")
load("sce.robj")
chord(seu = seu,sce=sce,doubletrat="estimated doubletrate 2",overkill=T,outname="the name you want 2",addmethods1 ="real_data.scores.csv",addmethods2 = "simulated_data.scores.csv" )
```

4.The doublets' barcodes are in the file "outname2_doublets.csv" 

## References

Allaire, J.J., et al. (2018). Reticulate: interface to Python. R Package Version 1.https://rstudio.github.io/reticulate/.

Al'Aref, S.J., Singh, G., Choi, J.W., Xu, Z., Maliakal, G., van Rosendael, A.R., Lee, B.C., Fatima, Z., Andreini, D., Bax, J.J., et al. (2020). A Boosted Ensemble Algorithm for Determination of Plaque Stability in High-Risk Patients on Coronary CTA. JACC Cardiovasc Imaging 13, 2162-2173.

Amezquita, R.A., Lun, A.T.L., Becht, E., Carey, V.J., Carpp, L.N., Geistlinger, L., Marini, F., Rue-Albrecht, K., Risso, D., Soneson, C., et al. (2020). Orchestrating single-cell analysis with Bioconductor. Nat Methods 17, 137-145.

Bais, A.S., and Kostka, D. (2020). scds: computational annotation of doublets in single-cell RNA sequencing data. Bioinformatics 36, 1150-1158.

Bernstein, N.J., Fong, N.L., Lam, I., Roy, M.A., Hendrickson, D.G., and Kelley, D.R. (2020). Solo: Doublet Identification in Single-Cell RNA-Seq via Semi-Supervised Deep Learning. Cell Systems 11, 95-101.e105.

DePasquale, E.A.K., Schnell, D.J., Van Camp, P.J., Valiente-Alandi, I., Blaxall, B.C., Grimes, H.L., Singh, H., and Salomonis, N. (2019). DoubletDecon: Deconvoluting Doublets from Single-Cell RNA-Sequencing Data. Cell Rep 29, 1718-1727 e1718.

Dietterich, T.G. (2000). Ensemble methods in machine learning. Paper presented at: International workshop on multiple classifier systems (Springer).

Fang, L.T., Afshar, P.T., Chhibber, A., Mohiyuddin, M., Fan, Y., Mu, J.C., Gibeling, G., Barr, S., Asadi, N.B., Gerstein, M.B., et al. (2015). An ensemble approach to accurately detect somatic mutations using SomaticSeq. Genome Biol 16, 197.

Fay, M.P., and Proschan, M.A. (2010). Wilcoxon-Mann-Whitney or t-test? On assumptions for hypothesis tests and multiple interpretations of decision rules. Stat Surv 4, 1-39.

Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A.K., Slichter, C.K., Miller, H.W., McElrath, M.J., Prlic, M., et al. (2015). MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biol 16, 278.

Grau, J., Grosse, I., and Keilwagen, J. (2015). PRROC: computing and visualizing precision-recall and receiver operating characteristic curves in R. Bioinformatics 31, 2595-2597.

Kang, H.M., Subramaniam, M., Targ, S., Nguyen, M., Maliskova, L., McCarthy, E., Wan, E., Wong, S., Byrnes, L., Lanata, C.M., et al. (2018). Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nat Biotechnol 36, 89-94.

Lambrechts, D., Wauters, E., Boeckx, B., Aibar, S., Nittner, D., Burton, O., Bassez, A., Decaluwé, H., Pircher, A., Van den Eynde, K., et al. (2018). Phenotype molding of stromal cells in the lung tumor microenvironment. Nat Med 24, 1277-1289.

Li, C., Liu, B., Kang, B., Liu, Z., Liu, Y., Chen, C., Ren, X., and Zhang, Z. (2020). SciBet as a portable and fast single cell type identifier. Nature communications 11, 1818.

Li, W.V., and Li, J.J. (2019). A statistical simulator scDesign for rational scRNA-seq experimental design. Bioinformatics 35, i41-i50.

Liu, B., Li, C., Li, Z., Wang, D., Ren, X., and Zhang, Z. (2020). An entropy-based metric for assessing the purity of single cell populations. Nature communications 11, 3155.

Lun, A.T., McCarthy, D.J., and Marioni, J.C. (2016). A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Res 5, 2122.

McGinnis, C.S., Murrow, L.M., and Gartner, Z.J. (2019). DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Systems 8, 329-337.e324.

McGinnis, C.S., Patterson, D.M., Winkler, J., Conrad, D.N., Hein, M.Y., Srivastava, V., Hu, J.L., Murrow, L.M., Weissman, J.S., Werb, Z., et al. (2019). MULTI-seq: sample multiplexing for single-cell RNA sequencing using lipid-tagged indices. Nat Methods 16, 619-626.

Potter, S.S. (2018). Single-cell RNA sequencing for the study of development, physiology and disease. Nat Rev Nephrol 14, 479-492.

Prakadan, S.M., Shalek, A.K., and Weitz, D.A. (2017). Scaling by shrinking: empowering single-cell 'omics' with microfluidic devices. Nat Rev Genet 18, 345-361.

Qiu, X., Mao, Q., Tang, Y., Wang, L., Chawla, R., Pliner, H.A., and Trapnell, C. (2017). Reversed graph embedding resolves complex single-cell trajectories. Nat Methods 14, 979-982.

Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J.C., and Müller, M. (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 12, 77.

Stoeckius, M., Zheng, S., Houck-Loomis, B., Hao, S., Yeung, B.Z., Mauck, W.M., 3rd, Smibert, P., and Satija, R. (2018). Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. Genome Biol 19, 224.

Street, K., Risso, D., Fletcher, R.B., Das, D., Ngai, J., Yosef, N., Purdom, E., and Dudoit, S. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics 19, 477.

Wolock, S.L., Lopez, R., and Klein, A.M. (2019). Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Syst 8, 281-291.e289.

Wu, Y., and Zhang, K. (2020). Tools for the analysis of high-dimensional single-cell RNA sequencing data. Nat Rev Nephrol 16, 408-421.

Xi, N.M., and Li, J.J. (2020). Benchmarking Computational Doublet-Detection Methods for Single-Cell RNA Sequencing Data. Cell Syst.

Yang, S., Corbett, S.E., Koga, Y., Wang, Z., Johnson, W.E., Yajima, M., and Campbell, J.D. (2020). Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biol 21, 57.

Zheng, G.X., Terry, J.M., Belgrader, P., Ryvkin, P., Bent, Z.W., Wilson, R., Ziraldo, S.B., Wheeler, T.D., McDermott, G.P., Zhu, J., et al. (2017). Massively parallel digital transcriptional profiling of single cells. Nature communications 8, 14049.

