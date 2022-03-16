\name{filter_gene}
\alias{filter_gene}
\title{filter_gene}
\description{To remove random noise corrputed genes by OPFSC.}
\usage{
        filter_gene(data, nBins = 20, plot_option = 1)
}

\arguments{
 \item{data}{The log2 transform count data.}
 \item{nBins}{The number of bins for generalized linear model construction.}
 \item{plot_option}{Option to show gene scatter plot and the contructed gene filtering baseline.}
 }
 \value{	
 \item{filter_gene}{The filtered data.}
}
 