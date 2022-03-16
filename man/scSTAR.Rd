\name{scSTAR}
\alias{scSTAR}
\title{scSTAR}
\description{The main function to extract the cell kinetics between to two experimental conditions, for example control and case. }

\usage{
scSTAR(Ctr, Case, anchorCells)
scSTAR(Ctr, Case, anchorCells, PLScomp = 0, NCV = 5, 
                       VG_Threshold = 0.01, BE_option = 1)
        }

\arguments{
 \item{Ctr, Case}{The OGFSC filtered data from the control and the case group, respectively. Genes vs. cells. }

\item{anchorCells}{The paired anchor cells identified from both control and case groups using the ¡®findAnchors¡¯ function.}

\item{PLScomp}{The number of PLS components, 0 - the value is estimated using cross-validation method; otherwise used the given fixed value. Be default 5.}

\item{NCV}{The number of folds of cross-validation. Be default 5.}

\item{VG_Threshold}{The p-value cutoff for DE gene identification.}
\item{BE_option}{The option to remove batch effect, 1 - perform batch effect removal (default), 0 - do not.}
}


 \value{	
 \item{Case_kinetics}{The scSTAR processed data from the case group.}
 \item{Ctr_kinetics}{The scSTAR processed data from the control group.}
 \item{idx_DE_up}{The index of the identified DE genes which have higher expression in the case group.}
 \item{idx_DE_down}{ The index of the identified DE genes which have higher expression in the control group.}
 \item{Ctr_corrected, Case_corrected}{The data with batch effect corrected only.}
}
