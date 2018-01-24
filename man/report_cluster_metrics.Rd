\name{report_cluster_metrics}
\alias{report_cluster_metrics}
\alias{calculate_cluster_metrics}
\alias{calculate_cluster_metrics_by_cluster}
\alias{plot_cluster_metrics}
\alias{value}
\alias{size}
\alias{L1}
\alias{Overall.mean}
\alias{mean.prom}

\title{
Reports BEARscc metrics for clusters.
}
\description{
To quantitatively evaluate the results, three metrics are calculated
    from the noise consensus matrix: 'stability' is the average
    frequency with which cells within a cluster associate with each
    other across simulated replicates; 'promiscuity' measures the
    association frequency between cells within a cluster and those
    outside of it; and 'score' is the difference between 'stability'
    and 'promiscuity'. Importantly, 'score' reflects the overall "robustness"
    of a cluster to technical variance. These metrics may be calculated
    on cell or cluster-wise basis; here, they are calculated cluster-wise.
}
\usage{
report_cluster_metrics(cluster_labels, consensus_matrix,
    weighted_mean = FALSE, plot = FALSE, file = "Rplot")
}

\arguments{
\item{cluster_labels}{
Cluster labels for each cell across various cluster numbers
    and the original clustering.
}
\item{consensus_matrix}{
A noise consensus output by \code{compute_consensus()}
}
\item{weighted_mean}{
A flag indicating whether to weigh observed clusters evenly
    or scale them by the number of samples in the cluster.
}
\item{plot}{
A flag to determine whether to plot the boxplot of cluster metrics
    evaluated from the noise consensus with root \code{file}.
}
\item{file}{
A string indicating the root desired for the resulting plots of the function.
}
}

\value{
A melted \code{data.frame} of BEARscc metrics for each cluster:
    \tabular{cll}{
    [,1] \tab "Cluster.identity" \tab The number of the cluster within the
        respective clustering. \cr
    [,2] \tab "Cluster.size" \tab Number of samples in the cluster.\cr
    [,3] \tab "Metric" \tab Whether the metric is the BEARscc score,
        promiscuity, or stability.\cr
    [,4] \tab "Value" \tab Value of the relevant BEARscc metric for the
        cluster in a clustering.\cr
    [,5] \tab "Clustering" \tab The clustering pertinant to the cell-wise
        metrics described.\cr
    [,6] \tab "Singlet"\tab A binary output concerning whether the cluster
        consists of a single sample.\cr
    [,7] \tab "Clustering.Mean"\tab The average of the respective metric
        across cells of the clsuter.\cr
    }
}


\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}


\examples{
data(analysis_examples)

cluster_scores.df <- report_cluster_metrics(BEARscc_clusts.df, noise_consensus,
    plot=TRUE, file="example")
cluster_scores.df
}

\keyword{ list }
