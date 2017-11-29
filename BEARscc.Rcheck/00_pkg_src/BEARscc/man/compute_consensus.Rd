\name{compute_consensus}
\alias{compute_consensus}
\alias{create_cm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute consensus matrix.
}
\description{
Computes the consensus matrix using a data.frame of cluster labels across different BEARscc simulated technical replicates.
}
\usage{
compute_consensus(cluster_labels)
}

\arguments{
  \item{cluster_labels}{
  A data.frame of labels assigned to each sample (rownames) across various simulated technical replicates designed by BEARscc (colnames).
  }
}
\details{
We provide a visual and quantitative representation of the clustering variation on a cell-by-cell level by using cluster labels to compute the number of times any given pair of cells associates in the same cluster; this forms the ‘noise consensus matrix’. Each element of this matrix represents the fraction of simulated technical replicates in which two cells cluster together (the ‘association frequency’), after using a clustering method of the user's choice to generate a data.frame of clustering labels. This consensus matrix may be used to compute BEARscc metrics at both the cluster and cell level.
}
\value{
  When the number of samples are \emph{n}, then the noise consensus resulting from this function is an \emph{n}x\emph{n} matrix describing the fraction of simulated technical replicates in which each cell of the experiment associates with another cell.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}

\section{A brief description of subfunctions}{
  \code{compute_consensus} relies on the following subfunction to compute the noise consensus. This function obtains all of the necessary information form the options of  \code{compute_consensus}. \cr
  \itemize{
    \item \code{names=rownames(cluster_labels)}
    \item \code{create_cm(cluster_labels, names)}
 }
}

\seealso{
 \itemize{
  \item \code{cluster_consensus()}
  \item \code{report_cluster_metrics()}
  \item \code{report_cell_metrics()}
 }
}
\examples{
data("analysis_examples")

noise_consensus <- compute_consensus(clusters.df)
noise_consensus
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ error }% __ONLY ONE__ keyword per line
