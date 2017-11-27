\name{report_cell_metrics}
\alias{report_cell_metrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reports BEARscc metrics for cells.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
report_cell_metrics(cluster_labels, consensus_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (cluster_labels, consensus_matrix)
{
    clusters <- rownames(cluster_labels)
    cluster_labels$X1 <- rep("1", dim(cluster_labels)[1])
    cluster_labels <- data.table(cluster_labels)
    cell.scores <- list()
    cell.scores <- lapply(cluster_labels, calculate_cell_metrics,
        consensus_matrix = consensus_matrix, clusters = clusters)
    cell.scores <- data.table(melt(cell.scores, id.var = c("rn",
        "cell", "cluster.size", "metric")))
    cell.scores$L1 <- as.factor(gsub("^X", "", cell.scores$L1))
    cell.scores <- cell.scores[, -5, with = FALSE]
    colnames(cell.scores) <- c("Cluster identity", "Cell", "Cluster size",
        "Metric", "Value", "Clustering")
    cell.scores
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
