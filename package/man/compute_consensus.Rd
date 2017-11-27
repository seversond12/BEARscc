\name{compute_consensus}
\alias{compute_consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute consensus matrix.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute_consensus(cluster_labels, plot = FALSE, file = "Rplot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster_labels}{
%%     ~~Describe \code{cluster_labels} here~~
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
function (cluster_labels, plot = FALSE, file = "Rplot")
{
    orig.clusts <- cluster_labels[, length(colnames(cluster_labels))]
    cluster_labels <- cluster_labels[, -length(colnames(cluster_labels))]
    cm <- lapply(cluster_labels, create_cm, names = rownames(cluster_labels))
    cm <- Reduce("+", cm)/dim(cluster_labels)[2]
    cm
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
