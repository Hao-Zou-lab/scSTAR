# scKinetics
Single-cell RNA sequencing has revolutionised our ability of interrogating gene expression at single-cell resolution, which allows detailed study of tumorigenesis and consequences of treatment. A critical component in understanding the impact of treatment is to measure changes in mRNA levels of individual cells between different conditions. Such measurements are often challenging due to practical limitation as each cell can only be profiled within a specific condition. Here we present, scKinetics, a novel algorithm that mathematically models the cell kinetics of interest by generating a latent virtual projection of each actual cell between conditions. Using simulated and experimental data, we showed scKinetics can reliably extract weak biological signals and promote biological discovery. It enabled the discovery of a new prognosis-associated cellular subtype in Head and neck cancer and an unknown aging-related bifurcation pattern of immune cells. Furthermore, its clinical potential is demonstrated in an immunotherapy study where therapy-response prediction accuracy was improved. 

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
This step can reduce the interferences of random noise by removing noise corrupted genes

OGFSC_idx_1 = OGFSC(data1, nBins = 30, alpha=0.5, plot_option = 1)$OGFSC_idx      
OGFSC_idx_2 = OGFSC(data2, nBins = 30, alpha=0.5, plot_option = 1)$OGFSC_idx      
idx_OGFSC = intersect(OGFSC_idx_1, OGFSC_idx_2)      
Ctr_filtered = data1[idx_OGFSC,]    
Case_filtered = data2[idx_OGFSC,]      

## identify anchor cells using KNN
The number indicates the number of nearest cells. By default 3. Although the users are free to change such number, tha authors recommond to keep it. 
                                                                            
anchorCells = findAnchors(Ctr_filtered, Case_filtered, 3) 

## single-cell kinetics estimation
The function to estimate individual cell kinetics. When PLScomp = 0, the number of PLScomp can be automatically estimated. But it can also be manually setted. The users are suggested to explore the values from 2-5.
  
out = scKinetics(Ctr_filtered, Case_filtered, anchorCells, PLScomp = 4) 
Case_kinetics = out$Case_kinetics
idx_DE_up = out$idx_DE_up

## auto clustering
Although the processed case data (out$Case_kinetics) and control data (out$Ctr_kinetics) could be analyzed by most of existing methods, we also incorporate a k-means based clustering method for data interpretation. The markerGenes are the list of DE genes (e.g., idx_DE_up) associated with each data cluster (e.g., Case_kinetics). 

K = 3 
out2 = autoClustering(Case_kinetics, K, geneList[idx_OGFSC], idx_DE_up) 
clusterIdx = out2$clusterIdx
markerGenes = out2$markerGenes









