\name{report_cell_metrics}
\alias{report_cell_metrics}
\alias{calculate_cell_metrics}
\alias{calculate_cell_metrics_by_cluster}
\title{
Reports BEARscc metrics for cells.
}
\description{
 To quantitatively evaluate the results, three metrics are calculated from the noise consensus matrix: ‘stability’ is the average frequency with which cells within a cluster associate with each other across simulated replicates; ‘promiscuity’ measures the association frequency between cells within a cluster and those outside of it; and ‘score’ is the difference between ‘stability’ and ‘promiscuity’. Importantly, ‘score’ reflects the overall “robustness” of a cluster and its constitutive samples to technical variance. These metrics may be calculated on cell or cluster-wise basis; here, they are calculated cell-wise.
 }

\usage{
report_cell_metrics(cluster_labels, consensus_matrix)
}

\arguments{
  \item{cluster_labels}{
  Cluster labels for each cell across various cluster numbers and the original clustering.
}
  \item{consensus_matrix}{
 A noise consensus output by \code{compute_consensus()}
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

cell_scores.dt <- report_cell_metrics(BEARscc_clusts.df, noise_consensus)
cell_scores.dt
}

\keyword{ list }

