\name{report_cell_metrics}
\alias{report_cell_metrics}
\alias{calculate_cell_metrics}
\alias{calculate_cell_metrics_by_cluster}

\alias{.}
\alias{rn}

\title{
Reports BEARscc metrics for cells.
}
\description{
To quantitatively evaluate the results, three metrics are calculated from
    the noise consensus matrix: 'stability' is the average frequency with
    which cells within a cluster associate with each other across simulated
    replicates; 'promiscuity' measures the association frequency between
    cells within a cluster and those outside of it; and 'score' is the
    difference between 'stability' and 'promiscuity'. Importantly, 'score'
    reflects the overall "robustness" of a cluster and its constitutive
    samples to technical variance. These metrics may be calculated on
    cell or cluster-wise basis; here, they are calculated cell-wise.
}

\usage{
report_cell_metrics(cluster_labels, consensus_matrix)
}

\arguments{
\item{cluster_labels}{
Cluster labels for each cell across various cluster numbers and
    the original clustering.
}
\item{consensus_matrix}{
A noise consensus output by \code{compute_consensus()}
}
}

\value{
A melted \code{data.frame} of BEARscc metrics for each cell:
\tabular{cll}{
    [,1] \tab "Cluster.identity" \tab The number of the cluster within the
        respective clustering \cr
    [,2] \tab "Cell" \tab The identifier of the sample in question.\cr
    [,3] \tab "Cluster.size" \tab Number of samples in the cluster.\cr
    [,4] \tab "Metric" \tab Whether the metric is the BEARscc Score,
        Promiscuity, or Stability. \cr
    [,5] \tab "Value" \tab Value of the relevant BEARscc metric for the cell
        in a given clustering. \cr
    [,6] \tab "Clustering"\tab The clustering pertinant to the cell-wise
        metrics described. \cr
    }
}

\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}

\examples{
data(analysis_examples)

cell_scores.df <- report_cell_metrics(BEARscc_clusts.df, noise_consensus)
cell_scores.df
}

\keyword{ list }

