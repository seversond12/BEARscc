\name{cluster_consensus}
\alias{cluster_consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cluster the consensus matrix.
}
\description{
We have found it useful to identify the optimal number of clusters in terms of resiliance to noise by examining these metrics by cutting hierarchical clustering dendograms of the noise consensus and comparing the results to the original clustering labels. To do this create a vector containing each number of clusters one wishes to examine (the function automatically determines the results for the dataset as a single cluster) and then cluster the consensus with this function.
}
\usage{
cluster_consensus(consensus_matrix, cluster_num, method = "complete")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{consensus_matrix}{
%%     ~~Describe \code{consensus_matrix} here~~
}
  \item{cluster_num}{
%%     ~~Describe \code{cluster_num} here~~
}
  \item{method}{
hierarchical clustering method to be used
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
function (consensus_matrix, cluster_num, method = "complete")
{
    cclust <- hclust(as.dist(abs(1 - t(consensus_matrix))), method = method)
    clusts <- data.frame(cutree(cclust, k = cluster_num))
    clusts
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
