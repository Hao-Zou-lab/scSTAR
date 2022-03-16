scSTAR <- function(Ctr, Case, anchorCells, PLScomp = 0, NCV = 5, 
                       VG_Threshold = 0.01, BE_option = 1)
{
  # function scSTAR extract cell kinetics.
  # 
  # Input: Ctr and Case - the data matrix from control and case groups, with
  #                       genes by cells
  #        anchorCells - the index of anchor cells
  #        PLScomp - the number of PLS components, 0 - the value is estimated
  #                   using cross-validation method; otherwise used the given fixed value.
  #        NCV  -  number of bins
  #        VG_Threshold  -  the threshold to identify variable genes
  #        BE_option  - option to remove batch effect, 1 - perform batch effect removal, 0 - do not
  # Output: idx_DE_up, idx_DE_down - the index of up and down regulated genes
  #         Case_kinetics, Ctr_kinetics - the data containing biological signal only
  #         Ctr_corrected, Case_corrected - the control and case data with batch effect corrected
  # Dr J Hao & Dr X Zou SJTU
  ## PLS1 contruction for batch effect removal
  if (BE_option == 1)
  {
    #prep, NCV, PLScomp_def, minNC
    MODEL = PLSconstruct(t(Ctr[,anchorCells[,1]]), t(Case[,anchorCells[,2]]), 
                         prep= 'mc', NCV, PLScomp, 1)
    # batch effect removal from the case group
    buffer2 = t(Case)
    buffer2 = buffer2-matrix(1,dim(buffer2)[1],1)%*%MODEL$mu
    buffer2 = buffer2-buffer2%*%MODEL$XL%*%ginv(MODEL$XL,tol = 0)
    buffer2 = buffer2+matrix(1,dim(buffer2)[1],1)%*%MODEL$mu

    
    # batch effect removal from the control group
    buffer1 = t(Ctr)
    buffer1 = buffer1-matrix(1,dim(buffer1)[1],1)%*%MODEL$mu
    buffer1 = buffer1-buffer1%*%MODEL$XL%*%ginv(MODEL$XL,tol = 0)
    buffer1 = buffer1+matrix(1,dim(buffer1)[1],1)%*%MODEL$mu
    
    # output the data with batch effect removed
    Ctr_corrected = t(buffer1)
    Case_corrected = t(buffer2)
  } else {
    Ctr_corrected = Ctr
    Case_corrected = Case
  }
  # non-DE genes for Null model construction
  P = matrix(1, dim(Case_corrected)[1], 1)
  for (i in 1:dim(Case_corrected)[1])
  {
    P[i] = t.test(Ctr_corrected[i,], Case_corrected[i,],
                  var.equal = T)$p.value
  }
 
  I = order(P,decreasing = T)
  idx_nonDE = I[1:floor(length(P)/5)]
  
  ## cell kinetics estimation by PLS2
    MODEL2 = PLSconstruct(t(Ctr_corrected), t(Case_corrected), prep = 'mc', NCV, PLScomp, 2)
    # SS (Sum of Square) of PLS2
    temp = t(Case_corrected)
    temp = temp-matrix(1,dim(temp)[1],1)%*%MODEL2$mu
    temp = temp%*%MODEL2$XL%*%ginv(MODEL2$XL,tol = 0)
    temp = temp+matrix(1,dim(temp)[1],1)%*%MODEL2$mu
    Case_kinetics = t(temp)
    
    
    temp = t(Ctr_corrected)
    temp = temp-matrix(1,dim(temp)[1],1)%*%MODEL2$mu
    temp = temp%*%MODEL2$XL%*%ginv(MODEL2$XL,tol = 0)
    temp = temp+matrix(1,dim(temp)[1],1)%*%MODEL2$mu
    Ctr_kinetics = t(temp) 
    
    
    ## DE gene identification
    # up-regulated genes
    FC = rowMeans(Case_kinetics)-rowMeans(Ctr_kinetics) # fold change 
    SSY = rowSums(Case_kinetics^2)
    SSY_nonDE = SSY[idx_nonDE]
    N = length(SSY_nonDE)                 # Number of 'Experiments' In Data Set
    SSY_nonDEMean = mean(SSY_nonDE)   # Mean Of All Experiments At Each Value Of ¡®x¡¯
    SSY_nonDESEM = sd(SSY_nonDE)/sqrt(N)  # Compute ¡®Standard Error Of The Mean¡¯ Of All Experiments At Each Value Of ¡®x¡¯
    CI = qt((1-VG_Threshold/2), N-1)   
    SSY_nonDECI95 = SSY_nonDEMean+(SSY_nonDESEM*CI)   # CI of the null model
    idx_DE_up = intersect(which(SSY>=SSY_nonDECI95), which(FC>0))  # up-regulated genes are SS>95% CI and FC>0
    
    # down-regulated genes
    SSY = rowSums(Ctr_kinetics^2) 
    SSY_nonDE = SSY[idx_nonDE] 
    N = length(SSY_nonDE)                                                # Number of ¡®Experiments¡¯ In Data Set
    SSY_nonDEMean = mean(SSY_nonDE)                                     # Mean Of All Experiments At Each Value Of ¡®x¡¯
    SSY_nonDESEM = sd(SSY_nonDE)/sqrt(N)                               # Compute ¡®Standard Error Of The Mean¡¯ Of All Experiments At Each Value Of ¡®x¡¯
    CI = qt((1-VG_Threshold/2), N-1)    
    SSY_nonDECI95 = SSY_nonDEMean+(SSY_nonDESEM*CI)   # CI of the null model
    idx_DE_down = intersect(which(SSY>=SSY_nonDECI95), which(FC>0))       # down-regulated genes are SS>95% CI and FC<0
    
    out <- list(Case_kinetics = Case_kinetics, idx_DE_up =idx_DE_up,
                Ctr_corrected =  Ctr_corrected, Case_corrected = Case_corrected, 
                Ctr_kinetics = Ctr_kinetics, idx_DE_down = idx_DE_down)
    return (out)
  }


  PLSconstruct <- function(data1, data2, prep, NCV, PLScomp_def, minNC)
  {
    # function PLSconstruct constructs PLS model based on the given data: data1 and data2, 
    #         and parameter specifications.
    # 
    # Input: data1 and data2 - the data matrix from control and case groups with genes by cells.
    #        prep - preprocessing methods.
    #        NCV  - number of folds of cross-validation.
    #        PLScomp_def - the number of PLS components, 0 - the value is automatically estimated;
    #                   otherwise used the given value.
    #        minNC - the minimum number of PLS components.
    # Output: model - the PLS model containing all PLS parameters and outputs
    
    # reshape data for PLSDA model construction
    X = rbind(data1, data2)
    Y = matrix(0, dim(X)[1],2)
    Y[seq_len(dim(data1)[1]),1] <- 1
    Y[(dim(data1)[1]+1):dim(Y)[1],2] <- 1
    
    maxNC = 5
    
    # scaling the data
    tmp <-dim(X)
    nsx <- tmp[1]
    nvx <- tmp[2]
    tmp <- dim(Y)
    nsy <- tmp[1]
    nvy <- tmp[2]
    
    muX = t(colMeans(X))
    muY = t(colMeans(Y))
    
    sigmaX = t(apply(X,2,sd))
    sigmaY = t(apply(Y,2,sd))
    
    model <- list(NULL)
    
    if (prep == 'no')
    {
      model$preprocessing<-'no'
    } else if (prep == 'mc')
    {
      #disp('Meancentering')
      model$preprocessing<-'mc'
      X<-X-matrix(rep(muX,each = nsx),nsx)
      Y<-Y-matrix(rep(muY,each = nsy),nsy)
    } else if (prep == 'uv')
    {
      # disp('Univariance scaling')
      model$preprocessing<-'uv'
      X<-X-matrix(rep(muX,each = nsx),nsx)
      Y<-Y-matrix(rep(muY,each = nsy),nsy)
      X<-X/matrix(rep(sigmaX,each = nsx),nsx)
      Y<-Y/matrix(rep(sigmaY,each = nsy),nsy)
    } else if (prep == 'pa')
    {
      model$preprocessing<-'pa'
      X<-X-matrix(rep(muX,each = nsx),nsx)
      Y<-Y-matrix(rep(muY,each = nsy),nsy)
      X<-X/matrix(rep(sqrt(sigmaX),each = nsx),nsx)
      Y<-Y/matrix(rep(sqrt(sigmaY),each = nsy),nsy)
    } else
    {
      model<-NULL
      stop('Unknown Preprocessing\n')
    }
    
    
    # PLS model construction
    block_num = floor(dim(X)[1]/NCV)
    Q2_ori_1 = 0
    Q2_ori_2 = 0
    nsx = dim(X)[1]
    if (PLScomp_def == 0)  # estimate the optimal number of PLScomp
    {
      PLScomp = 1
    for (nc in 1:maxNC)
    {
      Q2Yhatcum = matrix(0,1,NCV)
      for (cv in 1:NCV)
      {
        # cross validation
        idx_test = NCV*(1:block_num-1)+cv
        idx_tr = 1:nsx
        idx_tr = idx_tr[-idx_test]
        X_test = X[idx_test,]
        Y_test = Y[idx_test,]
        X_tr = X[idx_tr,]
        Y_tr = Y[idx_tr,]
        
        # model construct
        plsrDT = plsr(Y_tr ~ X_tr,nc)
        
        BETA <- plsrDT$coefficients[,,nc, drop=FALSE]
        # cumulative coefficients Intercept = T
        dB <- dim(BETA)
        dB[1] <- dB[1] + 1
        dnB <- dimnames(BETA)
        dnB[[1]] <- c("(Intercept)", dnB[[1]])
        BInt <- array(dim = dB, dimnames = dnB)
        BInt[-1,,] <-  BETA
        
        for (i in seq(along = nc))
          BInt[1,,i] <- plsrDT$Ymeans - plsrDT$Xmeans %*% BETA[,,i]
        
        BETA <- BInt[,,1]
        
        # predicting on the model
        PredY = cbind(matrix(1,dim(X_test)[1],1),X_test)%*%BETA
        # calculating Q2
        SSY = sum(Y_test^2)
        Q2Yhatcum[cv] = 1-sum((PredY-Y_test)*(PredY-Y_test))/SSY
      }
      Q2 = mean(Q2Yhatcum[!is.na(Q2Yhatcum)])
      
      # check if the current Q2 value is maximum
      if (Q2_ori_2>Q2 & Q2_ori_2>Q2_ori_1 & Q2_ori_1>Q2)
      {
        PLScomp = nc - 2  # this is the critical number of PLS components
        break
      } else {
        Q2_ori_2 = Q2_ori_1
        Q2_ori_1 = Q2
      }
    }
    PLScomp = max(PLScomp, minNC)
    } else
    {
      PLScomp = PLScomp_def
    }
    
    # model construction with the estimated number of PLS components
    plsestcom = plsr(Y ~ X ,PLScomp)
    XL = loadings(plsestcom)
    
    BETA <- plsestcom$coefficients[,,PLScomp, drop=FALSE]
    # cumulative coefficients Intercept = T
    dB <- dim(BETA)
    dB[1] <- dB[1] + 1
    dnB <- dimnames(BETA)
    dnB[[1]] <- c("(Intercept)", dnB[[1]])
    BInt <- array(dim = dB, dimnames = dnB)
    BInt[-1,,] <-  BETA
    
    for (i in seq(along = PLScomp))
      BInt[1,,i] <- plsestcom$Ymeans - plsestcom$Xmeans %*% BETA[,,i]
    
    BETA <- BInt[,,1]
    
    
    model$XL = matrix(XL, dim(XL)[1],dim(XL)[2])
    model$BETA = BETA
    
    if (prep == 'no')
    {
    } else if (prep == 'mc')
    {
      model$mu = muX
    } else if (prep == 'uv')
    {
      model$mu = muX
      model$sigmaX = sigmaX
    } else if (prep == 'pa')
    {
      model$mu = muX
      model$sigmaX = sigmaX
    }
    
    return(model)
  }
  
  
  
  
  
  
  
  
  
  
  