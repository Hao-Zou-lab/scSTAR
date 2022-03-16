\name{findAnchors}
\alias{findAnchors}
\title{findAnchors}
\description{To find anchor cells for PLS1 model construction using the KNN method.}
\usage{
        findAnchors(Ctr, Case, k = 3, dim = 5)
}

\arguments{
 \item{Ctr}{The OGFSC filtered control data.}
 \item{Case}{The OGFSC filtered case data.}
 \item{k}{The number of neighbours for KNN. The default value is 3.}
 \item{dim}{The number of principle components used for data normalization. The default value is 5.}
 }
 \value{	
 \item{anchorCells}{The list of anchor cells, paired with control (column 1) vs. case (column 2).}
}
 