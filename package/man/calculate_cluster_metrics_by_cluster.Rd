\name{calculate_cluster_metrics_by_cluster}
\alias{calculate_cluster_metrics_by_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{report_cluster_metrics()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculate_cluster_metrics_by_cluster(clusters, cluster_labels, consensus_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{
%%     ~~Describe \code{clusters} here~~
}
  \item{cluster_labels}{
%%     ~~Describe \code{cluster_labels} here~~
}
  \item{consensus_matrix}{
%%     ~~Describe \code{consensus_matrix} here~~
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
function (clusters, cluster_labels, consensus_matrix)
{
    in.cluster <- as.character(names(cluster_labels[cluster_labels ==
        clusters]))
    out.cluster <- as.character(names(cluster_labels[!(cluster_labels ==
        clusters)]))
    if (length(in.cluster) > 1) {
        inside.stability <- (sum(consensus_matrix[in.cluster,
            in.cluster, drop = FALSE]) - length(in.cluster))/(length(in.cluster)^2 -
            length(in.cluster))
    }
    else {
        inside.stability <- 1
    }
    promiscuity.matrix <- as.matrix(consensus_matrix[in.cluster,
        out.cluster, drop = FALSE])
    if (length(in.cluster) >= length(out.cluster)) {
        outside.promiscuity <- mean(promiscuity.matrix)
    }
    else {
        prom.mean.dt <- data.table(data.frame(colMeans(promiscuity.matrix)),
            keep.rownames = TRUE)
        colnames(prom.mean.dt)[2] <- "mean.prom"
        prom.mean.dt <- prom.mean.dt[order(-mean.prom)]
        prom.mean.dt <- prom.mean.dt[1:length(in.cluster), ]
        outside.promiscuity <- mean(promiscuity.matrix[, as.character(prom.mean.dt$rn),
            drop = FALSE])
    }
    score <- inside.stability - outside.promiscuity
    metrics.dt <- data.table(data.frame(score = score, Promiscuity = outside.promiscuity,
        Stability = inside.stability, size = length(in.cluster)))
    metrics.dt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
