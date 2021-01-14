library(R.matlab)
library(OGFSC)
library(FNN)
library(MASS)
library(tsne)
# source('D:/Work/NB/paperSections/tidy scKinetics/scStratify_0.1.37/scStratify/R/autoClustering.R')
# source('D:/Work/NB/paperSections/tidy scKinetics/scStratify_0.1.37/scStratify/R/findAnchors.R')
# source('D:/Work/NB/paperSections/tidy scKinetics/scStratify_0.1.37/scStratify/R/scKinetics.R')
library(scKinetics)

setwd("E:/algorithms/scIdentify/tidy scKinetics/R package")

data = readMat('data1.mat')
data1 <- as.matrix((data$data1))
data = readMat('data2.mat')
data2 <- as.matrix((data$data2))
geneList <- as.matrix(unlist(data$geneList))

## gene filtering
OGFSC_idx_1 = OGFSC(data1, nBins = 30, alpha=0.5, plot_option = 1)$OGFSC_idx
OGFSC_idx_2 = OGFSC(data2, nBins = 30, alpha=0.5, plot_option = 1)$OGFSC_idx
idx_OGFSC = intersect(OGFSC_idx_1, OGFSC_idx_2) 
Ctr_filtered = data1[idx_OGFSC,] 
Case_filtered = data2[idx_OGFSC,] 

## identify anchor cells using KNN
anchorCells = findAnchors(Ctr_filtered, Case_filtered, 3) 

## single-cell kinetics estimation
out = scKinetics(Ctr_filtered, Case_filtered, anchorCells, PLScomp = 4) 
Case_kinetics = out$Case_kinetics
idx_DE_up = out$idx_DE_up
## auto clustering
K = 3 
out2 = autoClustering(Case_kinetics, K, geneList[idx_OGFSC], idx_DE_up) 
clusterIdx = out2$clusterIdx
markerGenes = out2$markerGenes
