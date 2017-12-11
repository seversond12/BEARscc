\name{analysis_examples}
\alias{analysis_examples}
\alias{noise_consensus}
\alias{cluster.list}
\alias{sim_replicates}
\alias{estimated_noise}
\alias{recluster}
\alias{clusters.df}
\alias{BEARscc_clusts.df}
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
at https://bitbucket.org/bsblabludwig/bearscc.git.
}
\usage{data("analysis_examples")}

\value{
A long list as described below in the format section.
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
    list \tab "estimated_noise" \tab Pre-processed result from running 
        \code{estimate_noiseparameters()} on the toy datasets found in the 
        data in \code{BEARscc_examples}.\cr
    \tab [1] \tab genewise: A data.frame listing gene-wise parameters 
        necessary for computing drop-out recover probabilities in order to 
        create the drop-out model. \cr
    \tab [2] \tab inferred_prob: A data.table of the estimated probabilities 
        for drop-out recovery for each gene and relevant value k.\cr
    \tab [3] \tab bayes_parameters: A data.table of the estimated drop-out 
        model paramters utilized by \code{create_noiseinjected_counts()}.\cr
    \tab [4] \tab ERCC_parameters: A data.frame of the estimated noise model 
        parameters utilized by \code{create_noiseinjected_counts()}.\cr
    \tab [5] \tab models.dt: A melted data.table of observed counts, 
        pure Poisson, pure Neg. Binomial, and the BEARscc models that allows 
        plotting to asses fit.\cr
    \tab [6] \tab original.counts: A data.frame of the original counts matrix 
        identical to \code{data.counts.df}.\cr
    \tab [7] \tab spike.conc: A data.frame listing actual transcript number 
        for each spike-in. \cr
    \tab [8] \tab spike.counts: A data.frame listing the observed spike-in 
        counts \cr
    list \tab "sim_replicates" \tab Toy list of previously computed output of 
        \code{create_noiseinjected_counts()}.\cr
    \tab [,1] \tab Counts data.frame of simulated replicate 1. \cr
    \tab [,2] \tab Counts data.frame of simulated replicate 2.\cr
    \tab [,3] \tab Counts data.frame of simulated replicate 3.\cr
    \tab [,4] \tab Counts data.frame of simulated replicate 4.\cr
    \tab [,5] \tab Counts data.frame of simulated replicate 5.\cr
    \tab [,6] \tab Counts data.frame of simulated replicate 6.\cr
    \tab [,7] \tab Counts data.frame of simulated replicate 7.\cr
    \tab [,8] \tab Counts data.frame of simulated replicate 8.\cr
    \tab [,9] \tab Counts data.frame of simulated replicate 9.\cr
    \tab [,10] \tab Counts data.frame of simulated replicate 10.\cr
    \tab [,11] \tab Counts data.frame of observed data.\cr
        matrix \tab "noise_consensus" \tab 50 by 50 matrix of previously 
        computed output from \code{compute_consensus()}, which is the 
        noise consensus matrix resulting from the toy dataset, 
        \code{BEARscc_examples}.\cr
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
