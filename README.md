# scKinetics
The R package source code and a demo script of scKinetics

# Step 1: intallation
1. install the OGFSC R package from https://github.com/XZouProjects/OGFSC-R, and all the associated packages.
2. Installing scKinetics from github:
dowaload the file 'scKinetics_0.1.3.tar.gz' and install the package from local path.

# Step 2: prepare testing data
download demo data: data1 and data2

# Step 3: run scKinetics
library(R.matlab)
library(OGFSC)
library(FNN)
library(MASS)
library(tsne)
library(scKinetics)

## load data
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
 
The demo script includes OGFSC based gene filtering, scKinetics data processing and k-means clustering on the processed data.
