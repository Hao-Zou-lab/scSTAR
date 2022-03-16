\name{autoClustering}
\alias{autoClustering}
\title{autoClustering}
\description{Clustering the scSTAR processed data using k-means and assign the DE genes to each cluster.}

\usage{
    autoClustering(data, K, geneList, index_DE)
}

\arguments{
\item{data}{The scSTAR processed data from one group. Genes vs. cells.}

\item{K}{The number of clusters.}

\item{geneList}{The list gene symble IDs in the data.}

\item{index_DE}{The index of the DE, which have higher expression in this group than the opposite group.}
}


\value{	

 \item{clusterIdx_kmeans}{The cluster labels of each cell.}

 \item{markerGenes}{The grouped DE genes assigned to each cell cluster. This output may be applied for clustering results interpretation.}
}



