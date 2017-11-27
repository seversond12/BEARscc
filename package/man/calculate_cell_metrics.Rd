\name{calculate_cell_metrics}
\alias{calculate_cell_metrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{report_cell_metrics()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculate_cell_metrics(cluster_labels, consensus_matrix, clusters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster_labels}{
%%     ~~Describe \code{cluster_labels} here~~
}
  \item{consensus_matrix}{
%%     ~~Describe \code{consensus_matrix} here~~
}
  \item{clusters}{
%%     ~~Describe \code{clusters} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cluster_labels, consensus_matrix, clusters)
{
    names(cluster_labels) <- clusters
    cluster_names <- as.vector(as.character(unique(cluster_labels)))
    cluster_names <- as.character(unique(cluster_labels))
    cell_metric_list <- data.table(melt(lapply(cluster_names,
        calculate_cell_metrics_by_cluster, consensus_matrix = consensus_matrix,
        cluster_labels = cluster_labels), id.vars = c("cell",
        "cluster.size")))
    colnames(cell_metric_list)[c(3, 5)] <- c("metric", "rn")
    cell_metric_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
