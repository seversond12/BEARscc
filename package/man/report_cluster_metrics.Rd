\name{report_cluster_metrics}
\alias{report_cluster_metrics}
\alias{calculate_cluster_metrics}
\alias{calculate_cluster_metrics_by_cluster}

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
  A melted data.table of BEARscc metrics for each cluster:
\tabular{cll}{
[,1] \tab "Cluster identity" \tab The number of the cluster within the respective clustering. \cr
[,2] \tab "Cluster size" \tab Number of samples in the cluster.\cr
[,3] \tab "Metric" \tab Whether the metric is the BEARscc score, promiscuity, or stability.\cr
[,4] \tab "Value" \tab Value of the relevant BEARscc metric for the cluster in a clustering.\cr
[,5] \tab "Clustering" \tab The clustering pertinant to the cell-wise metrics described.\cr
[,6] \tab "Singlet"\tab A binary output concerning whether the cluster consists of a single sample.\cr
[,7] \tab "Clustering Mean"\tab The average of the respective metric across cells of the clsuter.\cr
  }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
\packageAuthor{BEARscc}

Maintainer: \packageMaintainer{BEARscc}
}


\section{A brief description of subfunctions}{
  \code{report_cluster_metrics} relies on the following subfunctions to report BEARscc metrics for each cluster. These subfunctions share many common options with the user interactive function, \code{report_cluster_metrics}.For further detail please examine source code in the R directory of this package:\cr
  \itemize{
    \item \code{cluster_names<-rownames(cluster_labels)}
    \item \code{cluster_labels<-data.table(cluster_labels)}
    \item \code{cluster.scores<-list()}
    \item \code{oneclust.stability<-(sum(consensus_matrix)-dim(consensus_matrix)[1])/dim(consensus_matrix)[1]^2}
    \item \code{cluster.scores$X1<-melt(data.table(rn="1",score=oneclust.stability,Promiscuity=0,Stability=oneclust.stability, size=dim(consensus_matrix)[1]), id.vars=c("rn", "size"))}
    \item \code{colnames(cluster.scores$X1)[3]<-"metric"}
    \item \code{cluster.scores<-append(cluster.scores,lapply(cluster_labels, `calculate_cluster_metrics`,consensus_matrix, cluster_names))}
    \item \code{names(cluster_labels)<-cluster_names}
    \item \code{clusters<-as.vector(as.character(unique(cluster_labels)))}
    \item \code{names(clusters)<-as.character(unique(cluster_labels))}
    \item \code{metric_list<-lapply(clusters, `calculate_cluster_metrics_by_cluster`, consensus_matrix=consensus_matrix, cluster_labels=cluster_labels)}
  }
}

\examples{
data(analysis_examples)

cluster_scores.dt <- report_cluster_metrics(BEARscc_clusts.df,noise_consensus, plot=TRUE, file="example")
cluster_scores.dt
}

\keyword{ list }
