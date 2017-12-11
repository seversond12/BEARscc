\name{cluster_consensus}
\alias{cluster_consensus}


\title{
Cluster the consensus matrix.
}
\description{
This function will perform hierarchical clustering on the noise consensus matrix
    allowing the user to investigate the appropriate number of clusters,
    \emph{k}, considering the noise within the experiment.
}
\usage{
cluster_consensus(consensus_matrix, cluster_num, method = "complete")
}

\arguments{
\item{consensus_matrix}{
A noise consensus output by \code{compute_consensus()}.
}
\item{cluster_num}{
The number of clusters expected from the hierarchical clustering of
    the noise consensus matrix.
}
\item{method}{
The hierarchical clustering method to be used on the consensus.
}
}

\details{
We have found it useful to identify the optimal number of clusters in terms
    of resiliance to noise by examining these metrics by cutting hierarchical
    clustering dendograms of the noise consensus and comparing the results to
    the original clustering labels. To do this create a vector containing
    each number of clusters one wishes to examine (the function automatically
    determines the results for the dataset as a single cluster) and then
    cluster the consensus with this function.

Frequently one will want to assess multiple possible cluster number
    situations at once. In this case it is recommended that one use a
    \code{lapply} in conjunction with a vector of all biologically
    reasonable cluster numbers to fulfill the task of attempting to
    identify the optimal cluster number.
}
\value{
The output is a vector of cluster labels based on hierarchical clustering
    of the noise consensus. In the event that a vector is supplied for
    number of clusters in conjunction with \code{lapply}, then the output
    is a data.frame of the cluster labels for each of the various number
    of clusters deemed biologically reasonable by the user.
}

\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}

\seealso{
    \code{compute_consensus}
    \code{report_cluster_metrics}
    \code{report_cell_metrics}
}
\examples{
data(analysis_examples)

vector <- seq(from=2, to=5, by=1)
BEARscc_clusts.df <- cluster_consensus(noise_consensus, vector)
BEARscc_clusts.df

}

\keyword{ cluster }
\keyword{ optimize }
