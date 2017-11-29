\name{report_cluster_metrics}
\alias{report_cluster_metrics}
\alias{calculate_cluster_metrics}
\alias{calculate_cluster_metrics_by_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reports BEARscc metrics for clusters.
}
\description{
  To quantitatively evaluate the results, three metrics are calculated from the noise consensus matrix: ‘stability’ is the average frequency with which cells within a cluster associate with each other across simulated replicates; ‘promiscuity’ measures the association frequency between cells within a cluster and those outside of it; and ‘score’ is the difference between ‘stability’ and ‘promiscuity’. Importantly, ‘score’ reflects the overall “robustness” of a cluster to technical variance. These metrics may be calculated on cell or cluster-wise basis; here, they are calculated cluster-wise.
}
\usage{
report_cluster_metrics(cluster_labels, consensus_matrix, weighted_mean = FALSE, plot = FALSE, file = "Rplot")
}

\arguments{
  \item{cluster_labels}{
  Cluster labels for each cell across various cluster numbers and the original clustering.
}
  \item{consensus_matrix}{
 A noise consensus output by \code{compute_consensus()}
}
  \item{weighted_mean}{
A flag indicating whether to weigh observed clusters evenly or scale them by the number of samples in the cluster.
}
  \item{plot}{
  A flag to determine whether to plot the boxplot of cluster metrics evaluated from the noise consensus with root \code{file}.
}
  \item{file}{
  A string indiciating the root desired for the resulting plots of the function.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}


\examples{
data(analysis_examples)

cluster_scores.dt <- report_cluster_metrics(BEARscc_clusts.df,noise_consensus, plot=TRUE, file="example")
clusters.dt
}

\keyword{ list }
