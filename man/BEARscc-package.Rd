\name{BEARscc-package}
\alias{BEARscc-package}
\alias{BEARscc}
\docType{package}
\title{
\packageTitle{BEARscc}
}
\description{
\packageDescription{BEARscc}
}

\details{
Single-cell transcriptome sequencing data are subject to substantial 
    technical variation and batch effects that can confound the 
    classification of cellular sub-types. Unfortunately, current 
    clustering algorithms don't account for this uncertainty. To 
    address this shortcoming, we have developed a noise perturbation 
    algorithm called BEARscc that is designed to determine the extent 
    to which classifications by existing clustering algorithms are 
    robust to observed technical variation.

BEARscc makes use of ERCC spike-in measurements to model technical variance 
    as a function of gene expression and technical dropout effects on lowly 
    expressed genes. In our benchmarks, we found that BEARscc accurately 
    models read count fluctuations and drop-out effects across transcripts 
    with diverse expression levels. Applying our approach to publicly 
    available single-cell transcriptome data of mouse brain and intestine, 
    we have demonstrated that BEARscc identified cells that cluster 
    consistently, irrespective of technical variation. For more details, 
    see the manuscript that is now available on bioRxiv.
}
\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}
\references{
Source code and README: <https://bitbucket.org/bsblabludwig/bearscc/overview>
Associated preprint: <https://www.biorxiv.org/content/early/2017/06/05/118919>
}

\keyword{ cluster }
\keyword{ optimize }
\keyword{ models }
\keyword{ robust }