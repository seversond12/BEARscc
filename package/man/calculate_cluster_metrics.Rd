\name{calculate_cluster_metrics}
\alias{calculate_cluster_metrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{report_cluster_metrics()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculate_cluster_metrics(cluster_labels, consensus_matrix, cluster_names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster_labels}{
%%     ~~Describe \code{cluster_labels} here~~
}
  \item{consensus_matrix}{
%%     ~~Describe \code{consensus_matrix} here~~
}
  \item{cluster_names}{
%%     ~~Describe \code{cluster_names} here~~
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
function (cluster_labels, consensus_matrix, cluster_names)
{
    names(cluster_labels) <- cluster_names
    clusters <- as.vector(as.character(unique(cluster_labels)))
    names(clusters) <- as.character(unique(cluster_labels))
    metric_list <- lapply(clusters, calculate_cluster_metrics_by_cluster,
        consensus_matrix = consensus_matrix, cluster_labels = cluster_labels)
    metric_list <- melt(metric_list, id.vars = c("size"))
    colnames(metric_list)[c(2, 4)] <- c("metric", "rn")
    metric_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
