\name{calculate_cell_metrics_by_cluster}
\alias{calculate_cell_metrics_by_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subfunction in \code{report_cell_metrics()}.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculate_cell_metrics_by_cluster(cluster_names, cluster_labels, consensus_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster_names}{
%%     ~~Describe \code{cluster_names} here~~
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
function (cluster_names, cluster_labels, consensus_matrix)
{
    in.cluster <- as.character(names(cluster_labels[cluster_labels ==
        cluster_names]))
    out.cluster <- as.character(names(cluster_labels[!(cluster_labels ==
        cluster_names)]))
    if (length(in.cluster) < length(out.cluster)) {
        if (length(in.cluster) == 1) {
            inside.stability <- 1
            outside.promiscuity <- max(consensus_matrix[in.cluster,
                out.cluster])
            score <- 1 - outside.promiscuity
            cell_metric.dt <- data.frame(cell = as.character(in.cluster[1]),
                score = score, Promiscuity = outside.promiscuity,
                Stability = inside.stability, cluster.size = length(in.cluster))
            cell_metric.dt$cluster.size <- length(in.cluster)
        }
        else {
            stability.dt <- data.table(data.frame(consensus_matrix[in.cluster,
                in.cluster, drop = FALSE]), keep.rownames = TRUE)
            stability.dt <- stability.dt[, .(Stability = (rowSums(.SD) -
                1)/(length(in.cluster) - 1)), keyby = rn]
            promiscuity.dt <- data.table(data.frame(consensus_matrix[in.cluster,
                out.cluster, drop = FALSE]), keep.rownames = TRUE)
            promiscuity.dt <- promiscuity.dt[, .(Promiscuity = apply(.SD,
                1, function(x, y = length(in.cluster)) {
                  mean(sort(x, decreasing = TRUE)[1:y])
                })), keyby = rn]
            cell_metric.dt <- stability.dt[promiscuity.dt]
            colnames(cell_metric.dt)[1] <- "cell"
            cell_metric.dt$score <- cell_metric.dt$Stability -
                cell_metric.dt$Promiscuity
            cell_metric.dt$cluster.size <- length(in.cluster)
        }
    }
    else if (length(in.cluster) >= length(out.cluster)) {
        if (length(in.cluster) == dim(consensus_matrix)[1]) {
            stability.dt <- data.table(data.frame(consensus_matrix[in.cluster,
                in.cluster, drop = FALSE]), keep.rownames = TRUE)
            stability.dt <- stability.dt[, .(Stability = (rowSums(.SD) -
                1)/(length(in.cluster) - 1)), by = rn]
            stability.dt$Promiscuity <- 0
            colnames(stability.dt)[1] <- "cell"
            cell_metric.dt <- stability.dt
            cell_metric.dt$cluster.size <- length(in.cluster)
            cell_metric.dt$score <- cell_metric.dt$Stability -
                cell_metric.dt$Promiscuity
        }
        else {
            stability.dt <- data.table(data.frame(consensus_matrix[in.cluster,
                in.cluster, drop = FALSE]), keep.rownames = TRUE)
            stability.dt <- stability.dt[, .(Stability = (rowSums(.SD) -
                1)/(length(in.cluster) - 1)), keyby = rn]
            promiscuity.dt <- data.table(data.frame(consensus_matrix[in.cluster,
                out.cluster, drop = FALSE]), keep.rownames = TRUE)
            promiscuity.dt <- promiscuity.dt[, .(Promiscuity = rowMeans(.SD)),
                keyby = rn]
            cell_metric.dt <- stability.dt[promiscuity.dt]
            cell_metric.dt$score <- cell_metric.dt$Stability -
                cell_metric.dt$Promiscuity
            cell_metric.dt$cluster.size <- length(in.cluster)
            colnames(cell_metric.dt) <- c("cell", "Stability",
                "Promiscuity", "score", "cluster.size")
        }
    }
    cell_metric.dt
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
