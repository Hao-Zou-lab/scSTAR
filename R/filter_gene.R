filter_gene <- function(data, nBins = 60, minBinSize = 100, LR_p = 0.01,
                    TW_threshold = 0.0001, plot_option = 0, paral_option = 1, 
                    CV = 10, maxPLScomp = 30, paralSize = 4, 
                    scalingMethod = c('mc', 'pa'), sampelLabels = NULL)
{
  #filter_gene to perform optimized gene filtering for single-cell RNA-seq data
  #Xin Zou, Jie Hao
  #2021
  #nBins: number of bins
  #minBinSize: minimeansDatam bin size
  #LR_p: p-value threshold to identify valid linear regression model
  #TW_threshold: the threshold to determine the number of eigenvalues on TW distribution
  #plot_option: plot the outputs of filter_gene, by default 0.
  
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  
  alpha <- 0.5 ###alpha = 0.5
  CVnum <- CV
  
  ## flag to select which algorithm is to be applied
  if (is.null(sampelLabels))
  {
    superviseFlag = 0
  } else {
    superviseFlag = 1
    if (paral_option == 1)
    {  
      cl<-makeCluster(paralSize, type = "SOCK")
      registerDoSNOW(cl)
    }
  }
  
  ## Pre-processing
  dummyM = matrix(0,dim(data)[2],2)
  dummyM[which(sampelLabels==1),1] = 1
  dummyM[which(sampelLabels==2),2] = 1
  
  ## remove non-sense or dropout values
  data_ori = data
  meansData = rowMeans(data)
  varsData = apply(data, 1, var)
  cv2Data = varsData/meansData^2
  idx = c(which(is.nan(cv2Data)), which(cv2Data==0), which(meansData==0))
  II = 1:nrow(data)
  if (length(idx)>0)
  {
  cv2Data = cv2Data[-idx]
  meansData= meansData[-idx]
  varsData= varsData[-idx]
  data=data[-idx,]
  II=II[-idx]
  }
  ## Binning genes and construct regression model using MLM method
  log10mean = log10(meansData)
  minV = min(log10mean)
  maxV = max(log10mean)
  stepsize = (maxV-minV)/nBins
  boundaries = seq(minV, maxV, stepsize)
  numBins = length(boundaries)-2
  Bins = matrix(0, numBins,2)
  BinSize = rep(0,numBins)
  BinIdx = vector('list', numBins)
  for (i in 1:numBins)
  {
    Bins[i,1] = boundaries[i]
    Bins[i,2] = boundaries[i+2]
    BinIdx[[i]] = which(log10mean>=Bins[i,1] & log10mean<Bins[i,2])
    BinSize[i] = length(BinIdx[[i]])
  }
  IdxValidBin = which(BinSize>=minBinSize)
  
  B = matrix(0,2,length(IdxValidBin))
  P = rep(1,length(IdxValidBin))
  for (i in 1:length(IdxValidBin))
  {
    idx_temp = BinIdx[[IdxValidBin[i]]]
    mean_temp = meansData[idx_temp]
    cv2_temp = cv2Data[idx_temp]
    #X = cbind(matrix(1,length(mean_temp), 1), 1/mean_temp)
    X = 1/mean_temp
    
    lm.model <- lm(cv2_temp ~ 1+X)
    B[,i] = matrix(coefficients(lm.model),2,1)
    if (dim(summary(lm.model)$coefficients)[1]<2)
    {
      P[i] = 1
    } else {
      P[i] = summary(lm.model)$coefficients[2,4]
    }
  }
  
  idx3 = which(P>LR_p)
  M = length(P)/2
  tmp = which(idx3<M)
  if (length(tmp) == 0)
  {  idx4 = idx3[tmp]
  }else
  {
    III = max(which(idx3<M))
    idx4 = idx3[III]
  }
  tmp <- which(idx3>M)
  if (length(tmp) == 0)
  {
    idx5 = idx3[tmp]
  } else
  {
    III = min(tmp)
    idx5 = idx3[III]
  }
  
  if ((length(idx4) != 0) && (length(idx5)!=0))
  {
    IdxValidBin <- IdxValidBin[-c(1:idx4,idx5:length(IdxValidBin))]
    B <-B[,-c(1:idx4,idx5:dim(B)[2])]
  }else if ((length(idx4) != 0) && (length(idx5)==0))
  {
    IdxValidBin <- IdxValidBin[-(1:idx4)]
    B <- B[,-(1:idx4)]
  } else if ((length(idx4) == 0) && (length(idx5)!=0))
  {
    IdxValidBin <-IdxValidBin[-(idx5:length(IdxValidBin))]
    B<-B[,-(idx5:dim(B)[2])]
  }
  
  idx = which(c(IdxValidBin,max(IdxValidBin))-c(0,IdxValidBin)>1)
  if (length(idx) != 0)
  {
    idx_temp = 1:min(idx)-1
    idx_temp = which(idx_temp>=1)
    if (length(idx_temp)>0)
    {
      IdxValidBin <- IdxValidBin[-idx_temp]
      B <- B[,-idx_temp]
    }
  }
  lowBoundary = Bins[min(IdxValidBin), 1]
  Bins = Bins[IdxValidBin,]
  BinIdx = BinIdx[IdxValidBin]
  BinSize = BinSize[IdxValidBin]
  idx4 <- NULL
  for (i in 1:dim(Bins)[1])
  {
    idx4 = c(idx4, BinIdx[[i]])
  }
  
  idx4 = sort(unique(idx4)) # the list of genes have valid expression levels
  BinAssignment = matrix(0, dim(Bins)[1], length(idx4)) # assign each gene to one or two bins, as bins have overlap
  for (i in 1:dim(Bins)[1])
  {
    Locb = match(BinIdx[[i]], idx4)
    BinAssignment[i,Locb] = 1 # 1 means assignment to the bin, 0 means no assignment
  }
  cv2_hat = rep(0,length(idx4))
  for (i in 1:length(idx4)) # estimate average noise cv2Data for each gene
  {
    idx5 = which(BinAssignment[,i]==1)
    temp = 0
    for (j in 1:length(idx5)) # for each LR model
    {
      temp = temp+B[1,idx5[j]]+B[2,idx5[j]]/meansData[idx4[i]]
    }
    cv2_hat[i] = temp/j
  }
  
  idx_output = II[idx4]
  cv2_output = cv2Data[idx4]
  
  ## plot regression curve if plot_option = 1
  if (plot_option == 1)
  {
    #x11()
    plot(log10(meansData), log10(cv2Data), type = "n" ,xaxt="n",yaxt="n",
         xlab = expression(bold(mu)), ylab= expression(CV^{2}),cex.lab=1.2)
    points(log10(meansData), log10(cv2Data), pch = ".", col = "gray",cex =1.2)
    points(log10(meansData[idx4]), log10(cv2_hat), pch = ".",col = "red", cex =2)
    xtick = floor(min(log10(meansData))):ceiling(max(log10(meansData)))
    xticklabel = (10^xtick)
    axis(1,xtick,labels=xticklabel,lwd = 2, cex.axis=1.5,cex.lab=3)
    
    ytick = floor(min(log10(cv2Data))):ceiling(max(log10(cv2Data)))
    yticklabel = (10^ytick)
    axis(2,ytick,labels=yticklabel,lwd = 2, cex.axis=1.5,cex.lab=3)
    
  }
  
 
  {
    # Using thresholding curves by alpha = 0.5
    data = data_ori[idx_output,]
    df = dim(data)[2]-1
    N = rep(0,length(alpha))
    selectedGenesIdx = vector('list', length(alpha))
    idx = which(cv2_output>cv2_hat*qchisq(alpha,df)/df)
    OGFSC_idx = idx_output[idx]
    cv2_threshold = cv2_hat*qchisq(alpha,df)/df

  
  
  if (!is.null(sampelLabels))
  {
    if (paral_option == 1)
    {  
      stopCluster(cl)
    }
  }
  
  OGF <- list(OGFSC_idx = OGFSC_idx, idx_output = idx_output)
  warnRead<-options(warn = warnDef)
  return(OGF)
  }
}

