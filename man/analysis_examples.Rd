\name{analysis_examples}
\alias{analysis_examples}
\alias{noise_consensus}
\alias{recluster}
\alias{clusters.df}
\alias{BEARscc_clusts.df}
\alias{BEAR_analyzed.sce}
\alias{.Random.seed}

\docType{data}

\title{
BEARscc downstream example objects.
}
\description{
The \code{analysis_examples} Rdata object contains downstream data objects
for use in various help pages for dynamic execution resulting from running
tutorial in README and vignette on \code{BEARscc_examples}. The objects
are a result of applying BEARscc functions as described in the README found
at https://bitbucket.org/bsblabludwig/bearscc.git or the vignette that
accompanies this package.
}
\usage{data("analysis_examples")}

\value{
An R data file with a \code{function}, \code{matrix}, two \code{data.frame}
    objects and a \code{SingleCellExperiment} object.
}

\format{
An R data file with: \cr
\tabular{ccl}{
    data.frame \tab "BEARscc_clusts.df" \tab Toy dataframe of previously
        computed cluster labels for each cell across various cluster numbers
        and the original clustering.\cr
    \tab [,1] \tab 2 cluster case\cr
    \tab [,2] \tab 3 cluster case\cr
    \tab [,3] \tab 4 cluster case\cr
    \tab [,4] \tab 5 cluster case\cr
    \tab [,5] \tab Original cluster case\cr
    data.frame \tab "clusters.df" \tab Toy dataframe of previously computed
        hierarchical clustering of 10 BEARscc simulated technical replicates as
        described in README and original hierarchical clustering,
        see \code{recluster} below. \cr
    \tab [,1] \tab Cluster labels from simulated replicate 1. \cr
    \tab [,2] \tab Cluster labels from simulated replicate 2.\cr
    \tab [,3] \tab Cluster labels from simulated replicate 3.\cr
    \tab [,4] \tab Cluster labels from simulated replicate 4.\cr
    \tab [,5] \tab Cluster labels from simulated replicate 5.\cr
    \tab [,6] \tab Cluster labels from simulated replicate 6.\cr
    \tab [,7] \tab Cluster labels from simulated replicate 7.\cr
    \tab [,8] \tab Cluster labels from simulated replicate 8.\cr
    \tab [,9] \tab Cluster labels from simulated replicate 9.\cr
    \tab [,10] \tab Cluster labels from simulated replicate 10.\cr
    \tab [,11] \tab Original cluster case\cr
    function \tab "recluster" \tab A function used to quickly illustrate
        replicate-wise clustering and the resulting consensus aspect of
        BEARscc. It is utilized in the README.\cr
    matrix \tab "noise_consensus" \tab 50 by 50 matrix of previously
        computed output from \code{compute_consensus()}, which is the
        noise consensus matrix resulting from the toy dataset,
        \code{BEARscc_examples}.\cr
    SCEList \tab "BEAR_analyzed.sce" \tab A \code{SingleCellExperiment} object
        discussed extensively in the vignette and ReadME.
    }
}

\source{
These data are the result of running the README,
    https://bitbucket.org/bsblabludwig/bearscc.git, on a subset of
    observations obtained by Drs. Michael White and Richard Owen in
    the Xin Lu Lab. Samples were sequenced by the Wellcome Trust Center
    for Genomics, Oxford, UK. The original data used to generate these
    objects are available in full with GEO accession number, GSE95155.
}
\references{
Source code and README: <https://bitbucket.org/bsblabludwig/bearscc/overview>
Associated preprint: <https://www.biorxiv.org/content/early/2017/06/05/118919>
}

\examples{
data(analysis_examples)
}

\keyword{datasets}
