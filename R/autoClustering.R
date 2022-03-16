autoClustering <- function(data, K, geneList, idx_DE)
{
# function Clustering the data into K clusters by k-means
# 
# Input: data - the data matrix with genes by cells
#        K - the number of clusters
#        geneList - the list of gene IDs.
#        idx_DE - the index of DE genes
# Output: clusterIdx_kmeans - the cluster labels
#         markerGenes - the marker genes assigned to each cluster

## tSNE scatter plot color coded by cluster labels
Color = brewer.pal(8,'Dark2') 
X = t(data[idx_DE,])
X = X - matrix(rep(rowMeans(X),dim(X)[2]),,dim(X)[2])

PCA<- prcomp(X)
SCORE <- PCA$x
clusterIdx_kmeans = kmeans(SCORE[,1:5], K, nstart = 20)$cluster
#clusterIdx_kmeans = kmeans(SCORE(:,1:5), K, 'replicates', 20, 'distance', 'sqeuclidean') 
ydata = tsne(X, initial_config = SCORE[,1:5],
             initial_dims = 5,perplexity = 30,k = 5)# k = 5

for (j in 1:K)
{
    II = which(clusterIdx_kmeans == j)
    if (j == 1)
    {      
      plot(ydata[,1], ydata[,2], type = "n", pch = 19, lwd = 4, 
           col = Color[j],xlab="tSNE_1",ylab="tSNE_2",cex.lab=1.2)
    }
    points(ydata[II,1], ydata[II,2], pch = 19, cex = .5, col = Color[j])
}

# hold on;
# [n,c] = hist3([ydata(:,1), ydata(:,2)]);
# contour(c{1},c{2},n');
# legend({'C1', 'C2', 'C3'});


## assign each DE gene to a cluster 
markerCluster = rep(0, length(idx_DE))
for (i in 1:length(idx_DE))
{  
temp = data[idx_DE[i],] 
expr = rep(0,K) 
for (j in 1:K)
{
  expr[j] = mean(temp[which(clusterIdx_kmeans==j)]) 
}
markerCluster[i] = which.max(expr) 
}
IX1 = NULL
IX2 = NULL 

for (i in 1:K)
{
  IX1 = c(IX1, (idx_DE[which(markerCluster==i)]))
  IX2 = c(IX2,  (which(clusterIdx_kmeans==i)))
}

buffer = data[IX1,] 
buffer2 = buffer[,IX2] 
buffer2 = buffer2 - matrix(rep(rowMeans(buffer2),dim(buffer)[2]),,dim(buffer)[2]) 
# clustergram(buffer2,'Cluster', 3, 'ColumnPDist', 'correlation')

markerGenes = vector("list", K) 
for (j in 1:K)
{
  markerGenes[[j]] = geneList[idx_DE[which(markerCluster==j)]] 
  # if(is.na(markerGenes[[j]]))
  # {
  #   markerGenes[[j]] = NULL 
  # }
}
aC <- list(clusterIdx_kmeans = clusterIdx_kmeans, markerGenes = markerGenes)
return(aC)
}
