
library(scSTAR)
library(tsne)
library(pls)

setwd("G:/algorithms/scIdentify/tidy scSTAR/R package")

extdir<-system.file("extdata",package="scSTAR")

data=readMat(paste0(extdir,'/demo_data.mat'))

data1 <- as.matrix((data$data1))#control data
    
data2 <- as.matrix((data$data2))#case  data

geneList <- as.matrix(unlist(data$geneList))

## gene filtering
data <- cbind(data1, data2)
idx_OGFSC = filter_gene(data, nBins = 20, plot_option = 1)$OGFSC_idx
Ctr_filtered = data1[idx_OGFSC,] 
Case_filtered = data2[idx_OGFSC,] 

## identify anchor cells using KNN
anchorCells = findAnchors(Ctr_filtered, Case_filtered, 3) 

## single-cell kinetics estimation
out = scSTAR(Ctr_filtered, Case_filtered, anchorCells, PLScomp = 4) 
Case_kinetics = out$Case_kinetics                                                   
idx_DE_up = out$idx_DE_up   

## auto clustering
K = 3                                                         
out2 = autoClustering(Case_kinetics, K, geneList[idx_OGFSC], idx_DE_up)                                                 
clusterIdx = out2$clusterIdx                                                                  
markerGenes = out2$markerGenes
