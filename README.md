![图片](https://user-images.githubusercontent.com/17633478/142828826-fa297984-5359-44ee-89d1-be8b30eb3398.png)

# scSTAR
Single-cell RNA sequencing has revolutionised our ability of interrogating gene expression at single-cell resolution, which allows detailed study of tumorigenesis and consequences of treatment. A critical component in understanding the impact of treatment is to measure changes in mRNA levels of individual cells between different conditions. Such measurements are often challenging due to practical limitation as each cell can only be profiled within a specific condition. Here we present, scSTAR, a novel algorithm that mathematically models the cell kinetics of interest by generating a latent virtual projection of each actual cell between conditions. Using simulated and experimental data, we showed scSTAR can reliably extract weak biological signals and promote biological discovery. It enabled the discovery of a new prognosis-associated cellular subtype in Head and neck cancer and an unknown aging-related bifurcation pattern of immune cells. Furthermore, its clinical potential is demonstrated in an immunotherapy study where therapy-response prediction accuracy was improved. 

# Step 1: intallation

1. install 'pls', 'tsne', 'R.matlab' R package from Cran.
 
2. Installing scSTAR from github:
dowaload the file 'scSTAR_0.1.1.0.tar.gz' and install the package from local path.

# Step 2: run scSTAR

    library(R.matlab)                                                                    
    library(tsne)                                           
    library(scSTAR)         
    library(pls)

## load data

    extdir<-system.file("extdata",package="scSTAR")
    
    data=readMat(paste0(extdir,'/demo_data.mat'))
   
    data1 <- as.matrix((data$data1))  
        
    data2 <- as.matrix((data$data2))   
    
    geneList <- as.matrix(unlist(data$geneList))

## gene filtering 
This step can reduce the interferences of random noise by removing noise corrupted genes

    data <- cbind(data1, data2)
    
    idx_OGFSC = filter_gene(data, nBins = 20, plot_option = 1)$OGFSC_idx    
    
    Ctr_filtered = data1[idx_OGFSC,]      
    
    Case_filtered = data2[idx_OGFSC,]                 
    
## identify anchor cells using KNN
The number indicates the number of nearest cells. By default 3. Although the users are free to change such number, tha authors recommond to keep it. 
                                                                            
    anchorCells = findAnchors(Ctr_filtered, Case_filtered, 3) 

## single-cell kinetics estimation
The function to estimate individual cell kinetics. When PLScomp = 0, the number of PLScomp can be automatically estimated. But it can also be manually setted. The users are suggested to explore the values from 3-5.
  
    out = scSTAR(Ctr_filtered, Case_filtered, anchorCells, PLScomp = 4)   
    
    Case_kinetics = out$Case_kinetics         
    
    idx_DE_up = out$idx_DE_up                                 

## data illustration (optional)
Although the processed case data (out$Case_kinetics) and control data (out$Ctr_kinetics) could be analyzed by most of existing methods, we also incorporate a k-means based clustering method for data interpretation. The markerGenes are the list of DE genes (e.g., idx_DE_up) associated with each data cluster (e.g., Case_kinetics). 

    K = 3     
    
    out2 = autoClustering(Case_kinetics, K, geneList[idx_OGFSC], idx_DE_up)  
    
    clusterIdx = out2$clusterIdx     
    
    markerGenes = out2$markerGenes









