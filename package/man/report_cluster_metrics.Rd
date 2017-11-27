\name{report_cluster_metrics}
\alias{report_cluster_metrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reports BEARscc metrics for clusters.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
report_cluster_metrics(cluster_labels, consensus_matrix, weighted_mean = FALSE, plot = FALSE, file = "Rplot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster_labels}{
%%     ~~Describe \code{cluster_labels} here~~
}
  \item{consensus_matrix}{
%%     ~~Describe \code{consensus_matrix} here~~
}
  \item{weighted_mean}{
%%     ~~Describe \code{weighted_mean} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
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
function (cluster_labels, consensus_matrix, weighted_mean = FALSE,
    plot = FALSE, file = "Rplot")
{
    cluster_names <- rownames(cluster_labels)
    cluster_labels <- data.table(cluster_labels)
    cluster.scores <- list()
    oneclust.stability <- (sum(consensus_matrix) - dim(consensus_matrix)[1])/dim(consensus_matrix)[1]^2
    cluster.scores$X1 <- melt(data.table(rn = "1", score = oneclust.stability,
        Promiscuity = 0, Stability = oneclust.stability, size = dim(consensus_matrix)[1]),
        id.vars = c("rn", "size"))
    colnames(cluster.scores$X1)[3] <- "metric"
    cluster.scores <- append(cluster.scores, lapply(cluster_labels,
        calculate_cluster_metrics, consensus_matrix = consensus_matrix,
        cluster_names = cluster_names))
    cluster.scores.m <- data.table(melt(cluster.scores, id.var = c("rn",
        "size", "metric")))
    cluster.scores.m <- cluster.scores.m[, -4, with = FALSE]
    cluster.scores.m$L1 <- gsub("^X", "", cluster.scores.m$L1)
    cluster.scores.m$singlet <- cluster.scores.m$size == 1
    cluster.scores.m$singlet <- factor(gsub(FALSE, "Cell number > 1",
        gsub(TRUE, "Cell number = 1", cluster.scores.m$singlet)),
        levels = c("Cell number > 1", "Cell number = 1"))
    if (weighted_mean == TRUE) {
        mean.overall <- cluster.scores.m[, sum(value * size/as.numeric(cluster.scores.m[L1 ==
            "1", ]$size[1])), by = c("L1", "metric")]
    }
    else {
        mean.overall <- cluster.scores.m[, mean(value), by = c("L1",
            "metric")]
    }
    colnames(mean.overall)[3] <- "Overall.mean"
    setkeyv(cluster.scores.m, c("L1", "metric"))
    setkeyv(mean.overall, c("L1", "metric"))
    cluster.scores.m <- cluster.scores.m[mean.overall]
    numbers <- as.character(sort(unique(suppressWarnings(na.omit(as.numeric(gsub("^X",
        "", cluster.scores.m$L1)))))))
    other <- unique(cluster.scores.m$L1[!(gsub("X", "", cluster.scores.m$L1) \%in\%
        numbers)])
    level <- as.character(c(numbers, other))
    cluster.scores.m$L1 <- factor(cluster.scores.m$L1, levels = level)
    if (plot == TRUE) {
        g <- ggplot(data = cluster.scores.m, aes(x = L1, y = value)) +
            geom_boxplot(color = "grey20", fill = "grey85") +
            guides(fill = "none", color = "none") + xlab("\nCluster number") +
            ylab("Cluster metric\n") + facet_wrap(~metric, ncol = 1,
            nrow = 3, scales = "free") + ggtitle("Cluster technical consensus scores for various clustering\n") +
            labs(alpha = "Cluster size") + geom_point(aes(x = factor(L1),
            y = Overall.mean), fill = "red", size = 2.5, color = "black",
            pch = 23) + annotate("segment", x = -Inf, xend = Inf,
            y = -Inf, yend = -Inf) + annotate("segment", x = -Inf,
            xend = -Inf, y = -Inf, yend = Inf) + theme(axis.text.x = element_text(angle = 90,
            hjust = 1))
        savefile <- paste(file, "_cluster_scores.pdf", sep = "")
        ggsave(savefile, plot = g, device = "pdf", width = 8.5,
            height = 11, units = "in")
    }
    colnames(cluster.scores.m) <- c("Cluster identity", "Cluster size",
        "Metric", "Value", "Clustering", "Singlet", "Clustering Mean")
    cluster.scores.m
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
