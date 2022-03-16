findAnchors <- function(Ctr, Case, k = 3, dime = 5)
{
  # This function find anchor cells using KNN method
  # Input: Ctr and Case - the data matrix from control and case groups, with
  #                       genes by cells
  #        k - the number of nearest neighbours in KNN
  #        dime - the number of PCs for KNN
  # Output: anchorCells - N by 2 matrix containing indexs of anchor cells
  #                       of control (1) and case (2) groups
  
  
  ## normalize data
  X = t(cbind(Ctr, Case))
  
  PCA<- prcomp(X)
  SCORE <- PCA$x
  
  XS = SCORE[1:dim(Ctr)[2],1:dime] 
  YS = SCORE[(dim(Ctr)[2]+1):dim(SCORE)[1],1:dime] 
  for (i in 1:dim(XS)[1])
  { 
    XS[i,] = XS[i,]/norm(XS, type = "2") 
  }
  
  for (i in 1:dim(YS)[1])
  {    
    YS[i,] = YS[i,]/norm(YS, type = "2") 
  }
  
  ## KNN on euclidean distance, which is equivalent to cosine distance with normalized data
  Idx1 <- knn(as.data.frame(XS),as.data.frame(YS),cl=factor(1:dim(XS)[1]),k = k) 
  Idx1 <- attributes(Idx1)$nn.index
  Idx2 <- knn(as.data.frame(YS),as.data.frame(XS),cl=factor(1:dim(YS)[1]),k = k) 
  Idx2 <- attributes(Idx2)$nn.index
  ## identify anchor cell pairs, which should be the nearest neighbours to each other
  anchorCells = matrix(NaN,max(dim(Idx1)[1],dim(Idx2)[1]), 2) 
  cellPairBuffer = matrix(0, max(dim(Idx1)[1]+dim(Idx2)[1]),1) 
  m = 1 
  n = 1 
  for (i in 1:dim(Idx1)[1])
  {
    for (j in 1:dim(Idx1)[2])
    { 
      if (length(which(i==Idx2[Idx1[i,j],])) != 0)
      { 
        if ((length(which(i==cellPairBuffer)) == 0) && (length(which(Idx1[i,j]==cellPairBuffer)) == 0))
        { 
          cellPairBuffer[m] = i
          cellPairBuffer[m+1] = Idx1[i,j] 
          anchorCells[n,2] = i  # Case group
          anchorCells[n,1] = Idx1[i,j] # Control group
          n = n+1 
          m = m+2 
        }
        break 
      }
    }
  }
  idx = which(is.nan(anchorCells[,1])) 
  anchorCells = anchorCells[-idx,]
  return(anchorCells)
}

